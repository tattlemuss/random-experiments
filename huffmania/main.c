/*
 * Small chunk of C code to generate Huffman encoding of a small data
 * block, and decode it again.
 * 
 * Mainly here for me to learn about the small details of doing the encoding.
 * 
 * Design goals: keep small and simple, and don't allocate memory.
 * Also, split all the phases of encoding into separate functions and
 * structs so they can be bolted together in different ways.
 * 
 * Implemented blind after reading the Wikipedia article on Huffman! 
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>		/* for memset */

const char g_testData[] = "Now is the winter of our discount tents.\0";

/* Specifies the number of tokens to allow.
 * We could use 257 if we wanted to add a terminator token.
*/
#define TOKEN_SIZE		(256U)

/* Value to represent used/unset token ID */
#define INVALID_TOKEN	((size_t)-1)

/* The max theoretical bit count for an encoded table is

    1 + 2 + 3 + ... n - 1 where n is the number of used tokens.

    i.e. 1/2n(n+1) - 1

    This would be in the case where each non-leaf used up
    one leaf token.

    This is impossible for 256 tokens since there will be
    numeric overflow in the frequency table (need to use
    fibonacci to generate it.)

    "Max enc bits: 65791 bytes: 8223"
*/

#define MAX_ENCODING_TABLE_BITS		((TOKEN_SIZE+1U)*(TOKEN_SIZE)-1U)
#define MAX_ENCODING_TABLE_BYTES	((MAX_ENCODING_TABLE_BITS + 7U) / 8U)

/******************************************************************************/
/* Defines the max number of nodes allowed in the Huffman tree.
 * Each time we combine two nodes, we add one more, so the worst case
 * for 4 tokens is our starting 4 + (2 + 1)
 */
#define FREQ_TREE_MAX_NODES      (TOKEN_SIZE*2U-1U)

/* This stores all the relative frequencies we want to adhere to,
* and a tree of left/right nodes which can be used to decode the tree.
*/
typedef struct _FreqTree
{
    size_t	m_count[FREQ_TREE_MAX_NODES];
    size_t	m_left[FREQ_TREE_MAX_NODES];
    size_t	m_right[FREQ_TREE_MAX_NODES];
    size_t  m_top;                 /* Index which represents the head node of the tree */
} FreqTree;

/******************************************************************************/
/* Represents the frequency tree values encoded into a flat
* bitfield which we can use to encode the output data with */
typedef struct _EncodedTable
{
    uint8_t 	m_packedData[MAX_ENCODING_TABLE_BYTES];	/* The encoded bits, all packed together. */
    size_t		m_offset[TOKEN_SIZE];	            /* Offset of each character's bits into m_encTable (in bits) */
    size_t		m_bitcount[TOKEN_SIZE];	            /* number of encoded bits for this character */
    size_t		m_currentOffset;		            /* current write offset in BITS when adding data into m_encTable */
} EncodedTable;

/******************************************************************************/
/* Simple structure to allow bitwise read/write into a buffer.
* Currently has no overflow handling, or callbacks to shell
  out to malloc or file storage */
typedef struct _Stream
{
    uint8_t m_buffer[8192];
    size_t  m_offset;        /* current offset in BITS */
} Stream;

/******************************************************************************/
/*  Reset frequency tree to initial values */
void freq_tree_init(FreqTree* pFreqTree)
{
    uint32_t i = 0U;
    for (i = 0U; i < FREQ_TREE_MAX_NODES; ++i)
    {
        pFreqTree->m_count[i] = 0U;
        pFreqTree->m_left[i]  = INVALID_TOKEN;
        pFreqTree->m_right[i] = INVALID_TOKEN;
    }
    pFreqTree->m_top = FREQ_TREE_MAX_NODES;
}

/* Run through a buffer of tokens and accumulate frequencies of each token used.
 * Can be called multiple times to accumulate over time */
void freq_tree_accumulate(FreqTree* pFreqTree, const char* pData, size_t len)
{
    /* Increment frequency on each token read */
    uint32_t i = 0U;
    for (i = 0U; i < len; ++i)
        ++pFreqTree->m_count[*pData++];
}

/* Debug function to list frequencies used. */
void freq_tree_dump(const FreqTree* pFreqTree)
{
    uint32_t i = 0U;
    for (i = 0U; i < FREQ_TREE_MAX_NODES; ++i)
        if (pFreqTree->m_count[i] != 0U)
            printf("Frequency for %u is %u\n", i, pFreqTree->m_count[i]);
}

/******************************************************************************/
/* Reset all data in the encoded table */
void encoded_table_init(EncodedTable* pTable)
{
    memset(pTable->m_offset,   0, sizeof(pTable->m_offset));
    memset(pTable->m_bitcount, 0, sizeof(pTable->m_bitcount));
    pTable->m_currentOffset = 0U;
}

void set_encoded_bit(EncodedTable* pTable, size_t bit_offset, uint8_t value)
{
    uint8_t* pWritePtr = &pTable->m_packedData[bit_offset / 8U];
    uint8_t shift = 7U - (bit_offset % 8U);
    *pWritePtr |= value << shift;
}

uint8_t read_encoded_bit(const EncodedTable* pTable, size_t bit_offset)
{
    const uint8_t* pReadPtr = &pTable->m_packedData[bit_offset / 8U];
    uint8_t shift = 7U - (bit_offset % 8U);
    return ((*pReadPtr) >> shift) & 1U;
}

/* Pulls out the smallest huffman frequency. Returns the index. 
   Linear scan, i.e. slow.
   Can be sped up by presorting the frequencies and keeping 2 lists
   (see Wikipedia article) */
size_t get_smallest_huffman(const FreqTree* pFreqTree, size_t scan_size)
{
    size_t small_val = INVALID_TOKEN;
    size_t small_ind = INVALID_TOKEN;
    uint32_t i = 0U;
    //assert(scan_size <= FREQ_TREE_MAX_NODES);
    for (i = 0U; i < scan_size; ++i)
    {
        size_t val = pFreqTree->m_count[i];
        if ((val) && (val < small_val))
        {
            small_val = val;
            small_ind = i;
        }
    }
    return small_ind;
}

void generate_huffman(FreqTree* pFreqTree)
{
    size_t new_index = TOKEN_SIZE;
    while (1)
    {
        size_t ind1 = get_smallest_huffman(pFreqTree, new_index);
        if (ind1 == INVALID_TOKEN)
            break;

        size_t val1 = pFreqTree->m_count[ind1];
        pFreqTree->m_count[ind1] = 0U;		/* hide from the table -- not needed any more */

        size_t ind2 = get_smallest_huffman(pFreqTree, new_index);
        if (ind2 == INVALID_TOKEN)
            break;
        size_t val2 = pFreqTree->m_count[ind2];
        pFreqTree->m_count[ind2] = 0U;

        /* Make a new combined entry node */
        pFreqTree->m_count[new_index] = val1 + val2;
        pFreqTree->m_left[new_index] = ind1;		
        pFreqTree->m_right[new_index] = ind2;
        ++new_index;		
    }	
    /* the top of the tree will always be the last index added */
    pFreqTree->m_top = --new_index;
}

/* Structure kept on stack during freq_tree_traverse(), in order
 * to keep a history of the bytes to write out. */
typedef struct _BitChain
{
    const struct _BitChain	*m_pParent;
    uint8_t					m_value;
} BitChain;

/* Traverse the tree in order, depth-first, to create the encoded
 * bits for each Huffman token, and fill it into EncodedTable.
 * 'onlyCount' denotes whether to just count the bits used (for debug)
 * or write them out to EncodedTable. We could do this with a callback.
 */
size_t freq_tree_traverse(EncodedTable* pEnc, const FreqTree* pFreqTree, size_t top, 
                          const BitChain* pChain, size_t bitcount, uint8_t onlyCount)
{
    /* Print out if we are a leaf */
    if (pFreqTree->m_left[top] == INVALID_TOKEN)
    {
        /* This is a leaf -- do something! */
        /* Scan upwards for the huffman code (backwards) */
        size_t size = 0U;
        if (onlyCount)
        {
            /* Don't do anything here other than fill the bitcount of the token */
            pEnc->m_bitcount[top] = bitcount;
        }
        else
        {
            printf("Encoding token: %u\n", top);
            /* Phase 2: actually write the bits */
            pEnc->m_offset[top] = pEnc->m_currentOffset;

            /* We need to write backwards, so start at the end */
            pEnc->m_currentOffset += bitcount;
            size_t write_offset = pEnc->m_currentOffset;
            while (pChain)
            {
                --write_offset;
                set_encoded_bit(pEnc, write_offset, pChain->m_value);
                pChain = pChain->m_pParent;
            }
        }
        return bitcount;
    }
    BitChain bc;
    bc.m_pParent = pChain;
    bc.m_value = 0U;
    size_t totalsize =  freq_tree_traverse(pEnc, pFreqTree, pFreqTree->m_left[top],  &bc, bitcount + 1U, onlyCount);
    bc.m_value = 1U;
    totalsize        += freq_tree_traverse(pEnc, pFreqTree, pFreqTree->m_right[top], &bc, bitcount + 1U, onlyCount);
    return totalsize;
}

/* Build up the encoding table from the current frequency tree. */
size_t create_encodings(EncodedTable* pEnc, const FreqTree* pFreqTree)
{    encoded_table_init(pEnc);

    /* Phase 1: fill out sizes and calc total size */
    /* This isn't needed any more, but useful as a check */
    size_t total_bitcount = freq_tree_traverse(pEnc, pFreqTree, pFreqTree->m_top, NULL, 0U, 0U);
    printf("Total bitcount of all codes: %u bits\n", total_bitcount);

    /* Phase 2: fill out encoded data */	
    freq_tree_traverse(pEnc, pFreqTree, pFreqTree->m_top, NULL, 0U, 1U);
}

/* Write out the encodings for tokens that were used */
size_t encoded_table_dump(const EncodedTable* pEnc)
{
    uint32_t i = 0U;
    for (i = 0U; i < TOKEN_SIZE; ++i)
    {
        if (pEnc->m_bitcount[i])
        {
            uint32_t startBit = pEnc->m_offset[i];
            uint32_t endBit   = startBit + pEnc->m_bitcount[i];
            printf("Token: %u -> ", i);
            while (startBit != endBit)
            {
                printf("%u", read_encoded_bit(pEnc, startBit));
                ++startBit;
            }				
            printf("\n");			
        }
    }
}

/* Reset the current pointer in the stream back to the start. */
void stream_reset(Stream* pStream)
{
    pStream->m_offset = 0U;
}

/* Put a bit in the stream and advance the write ptr */
void stream_write_bit(Stream* pStream, uint8_t value)
{
    uint8_t* pWritePtr = &pStream->m_buffer[pStream->m_offset / 8U];
    uint8_t shift = 7U - (pStream->m_offset % 8U);
    *pWritePtr |= value << shift;
    ++(pStream->m_offset);
}

/* Get current bit from the stream and advance the read ptr */
uint8_t stream_read_bit(Stream* pStream)
{
    uint8_t* pReadPtr = &pStream->m_buffer[pStream->m_offset / 8U];
    uint8_t shift = 7U - (pStream->m_offset % 8U);
    ++(pStream->m_offset);
    return ((*pReadPtr) >> shift) & 0x1;
}

/* Encode the data buffer into the stream, using the encoding data 
 * provided. */
void stream_encode(Stream* pStream, const char* pData, size_t len, const EncodedTable* pEnc)
{
    uint32_t i = 0U;
    for (i = 0U; i < len; ++i)
    {
        uint8_t org = *pData++;

        /* Look up in encoded table */
        uint32_t startBit = pEnc->m_offset[org];
        uint32_t endBit   = startBit + pEnc->m_bitcount[org];
        while (startBit != endBit)
        {
            stream_write_bit(pStream, read_encoded_bit(pEnc, startBit));
            ++startBit;
        }               
    }
}

/* Decode the stream data and print out the contents. */
void stream_decode(Stream* pStream, const FreqTree* pFreqTree)
{
    while (1)
    {
        uint32_t current = pFreqTree->m_top;
        do
        {
            uint8_t bit = stream_read_bit(pStream);
            if (bit)
                current = pFreqTree->m_right[current];
            else
                current = pFreqTree->m_left[current];
            // Keep going until at a leaf
        }
        while (pFreqTree->m_left[current] != INVALID_TOKEN);
        printf("%c ", current);
        if (!current)
            // Terminator
            break;
    }
    printf("\n");
}

int main( char** argv )
{
    FreqTree  tree;
    freq_tree_init(&tree);

    /* Scan buffer once to build up frequency info */
    freq_tree_accumulate(&tree, g_testData, sizeof(g_testData));
    freq_tree_dump(&tree);
    
    /* Make the tree from the current set of frequencies */
    generate_huffman(&tree);

    /* Now "flatten" the tree out into a set of bits for each used token */
    EncodedTable enc;
    create_encodings(&enc, &tree);
    encoded_table_dump(&enc);
    
    printf("Original stream size: %u bytes\n", sizeof(g_testData));

    /* OK, let's use the encoding table to pack some Huffman data */
    Stream encodedStream;
    stream_reset(&encodedStream);
    stream_encode(&encodedStream, g_testData, sizeof(g_testData), &enc);       
    printf("Encoded stream size: %u bits, %u bytes\n", encodedStream.m_offset, (encodedStream.m_offset + 7U) / 8U);

    printf("Decoding to stdout:\n");
    /* Go back to start of stream */
    stream_reset(&encodedStream);
    
    /* We need the original tree to decode, scanning down into its leaves.
     * In the real world we'd need to serialise the tree
     */
    stream_decode(&encodedStream, &tree);
    printf("Decoded!\n");
    
    return 0;
}
