//g++ -o hirgc hirgc.cpp -std=c++17
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <cstdint>
#include <chrono>
using namespace std;

int K = 32;
int MAX_CANDIDATES = 32;

uint64_t encodeKmer(const string& str, size_t start) {
    uint64_t val = 0; // A i N se spremaju u 0
    for (int i = 0; i < K; i++) {
        val <<= 2;
        char c = str[start + i];

        if      (c == 'C') val |= 1;
        else if (c == 'G') val |= 2;
        else if (c == 'T') val |= 3;
        // inace           val |= 0;
    }
    return val;
}

struct Block {
    size_t pos;
    size_t len;
};

string readGenomeFile(string filename, string& metadata, vector<Block>* lowerBlocks = nullptr)
{
    ifstream f(filename);

    if (!f.is_open()) {
        cout << "ERROR: could not open file: " << filename << endl;
        exit(1);
    }

    string result = "";   
    string line = "";

    //zapise metadata s prvog retka .fa datoteke
    getline(f, metadata);

    while (getline(f, line))
    {
        if (line.length() == 0 || line[0] == '>')
            continue;

        for (int i = 0; i < (int)line.length(); i++) {
            if (lowerBlocks != nullptr && islower(line[i])) {
                if (!lowerBlocks->empty() && lowerBlocks->back().pos + lowerBlocks->back().len == result.length() + i) {
                    lowerBlocks->back().len++;
                } else {
                    lowerBlocks->push_back({result.length() + i, 1});
                }
            }
            line[i] = toupper(line[i]);
        }
        result += line; 
    }

    f.close();
    return result; // result = AGCTACGTACGT..
}

struct FlatIndex {
    uint32_t table_size;
    uint32_t table_mask;
    vector<uint32_t> head;
    vector<uint32_t> next;
};

FlatIndex buildIndex(const string& ref)
{   // Koristimo "parallel arrays" pristup
    FlatIndex index;
    index.table_size = 1 << 26; // 67 milijuna brzih slotova
    index.table_mask = index.table_size - 1;

    if (ref.length() < K) return index;

    // Alokacija dvije velike strukture u memoriji odjednom
    index.head.assign(index.table_size, (uint32_t)-1);
    index.next.assign(ref.length(), (uint32_t)-1);

    uint64_t mask = (K == 32) ? ~0ULL : ((1ULL << (K * 2)) - 1); // npr. 1000000 - 1 = 0111111 za K=3
    uint64_t kmer = encodeKmer(ref, 0) & mask;
    
    uint64_t hash_val = (kmer * 0xbf58476d1ce4e5b9ULL) & index.table_mask; // isto kao % table_size
    index.next[0] = index.head[hash_val];
    index.head[hash_val] = 0;

    size_t stopAt = ref.length() - K;

    // svaki k-mer se enkodira i sprema
    for (size_t i = 1; i <= stopAt; i++)
    {
        kmer <<= 2;
        char c = ref[i + K - 1];
        if (c == 'C') kmer |= 1;
        else if (c == 'G') kmer |= 2;
        else if (c == 'T') kmer |= 3;
        kmer &= mask;

        hash_val = (kmer * 0xbf58476d1ce4e5b9ULL) & index.table_mask; // 
        index.next[i] = index.head[hash_val];
        index.head[hash_val] = i;
    }

    cout << "index built." << endl;

    return index;
}

void compressGenome(string refFile, string targetFile, string outputFile)
{
    string refMetadata;
    cout << "reading reference genome..." << endl;
    string ref = readGenomeFile(refFile, refMetadata); // ref = "ACGTACGATATCATCG..."
    cout << "reference length: " << ref.length() << " bases" << endl;

    string targetMetadata;
    cout << "reading target genome..." << endl;
    vector<Block> lowerBlocks;
    string target = readGenomeFile(targetFile, targetMetadata, &lowerBlocks); // target = "CAGTCAGCAGGACGT..."
    if (targetMetadata.empty()) {
        targetMetadata = ">recovered_sequence";
    }
    cout << "target length: " << target.length() << " bases" << endl;
    cout << "building kmer index with K=" << K << "..." << endl;


    FlatIndex index = buildIndex(ref);

    // za zapamtit gdje su bili N znakovi
    vector<Block> nBlocks;
    
    size_t i = 0;
    while (i < target.length()) {
        if (target[i] == 'N') {
            size_t start = i;
            while (i < target.length() && target[i] == 'N') {
                target[i] = 'A'; // preimenujemo N u A da zapis ostane 2-bitan, vracamo u dekompresiji
                i++;
            }
            nBlocks.push_back({start, i - start});
        } else {
            i++;
        }
    }

    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cout << "ERROR: cant open output file " << outputFile << endl;
        exit(1);
    }


    outFile << "HIRG" << "\n";
    outFile << targetMetadata << "\n"; //metadata
    outFile << refFile << "\n"; 
    outFile << target.length() << "\n"; //duljina
    outFile << K << "\n";           

    outFile << nBlocks.size() << "\n"; // zapisuje gdje su bili N znakovi
    for (const auto& block : nBlocks) {
        outFile << block.pos << " " << block.len << "\n";
    }

    outFile << lowerBlocks.size() << "\n"; // zapisuje gdje su bila lower case slova
    for (const auto& block : lowerBlocks) {
        outFile << block.pos << " " << block.len << "\n";
    }

    size_t pos = 0;
    size_t numMatches  = 0;
    size_t numLiterals = 0;
    size_t totalMatchLen = 0;  



    while (pos < target.length())
    {
        bool canLookup = (pos + K) <= target.length();

        if (canLookup)
        {

            uint64_t kmer = encodeKmer(target, pos); //pritvori u dekadski broj
            uint64_t hash_val = (kmer * 0xbf58476d1ce4e5b9ULL) & index.table_mask;

            uint32_t bestRefPos = 0;
            size_t bestMatchLen = 0;
            
            uint32_t candidate = index.head[hash_val];
            int count = 0;

            while (candidate != (uint32_t)-1 && count < MAX_CANDIDATES)
            {
                //trazimo najdulji match
                size_t matchLen = 0;

                while (matchLen < 65535
                    && (pos + matchLen) < target.length()
                    && (candidate + matchLen) < ref.length()
                    && target[pos + matchLen] == ref[candidate + matchLen])
                {
                    matchLen++;
                }

                if (matchLen >= K && matchLen > bestMatchLen)
                {
                    bestMatchLen = matchLen;
                    bestRefPos   = candidate;
                }
                
                candidate = index.next[candidate];
                count++;

            }

            if (bestMatchLen >= K) //ako je match pronađen
            {
                // sprema "M index_pozicije duljina"
                outFile << "M " << bestRefPos << " " << bestMatchLen << "\n";

                pos += bestMatchLen;

                numMatches++;
                totalMatchLen += bestMatchLen;

                continue;
            }
        }
        
        // ako nije pronasao nijedan match u index-u, zapise literal
        outFile << "L " << target[pos] << "\n";
        pos++;      
        numLiterals++;
    }

    outFile.close();

    cout << "\n=== primary compression done ===" << endl;
    cout << "match records:    " << numMatches << endl;
    cout << "literal records:  " << numLiterals << endl;
    cout << "bases via matches:" << totalMatchLen << " / " << target.length() << endl;
    cout << "\nApplying 7-Zip compression to the output file..." << endl;
    
    // 7-zip kompresija
    string cmd = "7z a " + outputFile + ".7z " + outputFile;
    int result = system(cmd.c_str());

    string finalFile = outputFile;
    if (result == 0) {
        cout << "7-Zip compression successful. Created: " << outputFile << ".7z" << endl;
        finalFile = outputFile + ".7z";

    } else {
        cout << "Warning: 7-Zip compression failed." << endl;
    }



    ifstream checkSize(finalFile, ios::binary | ios::ate);
    long long compressedSize = checkSize.is_open() ? (long long)checkSize.tellg() : 0LL;
    if (checkSize.is_open()) checkSize.close();

    cout << "\n=== final stats ===" << endl;
    cout << "original size:    " << target.length() << " bytes" << endl;
    cout << "compressed size:  " << compressedSize << " bytes (" << finalFile << ")" << endl;

    if (target.length() > 0) {
        double ratio = (double)compressedSize / (double)target.length() * 100.0;
        cout << "compression ratio: " << ratio << "%" << endl;
    }
}

int main(int argc, char* argv[])
{
    auto start_time = chrono::high_resolution_clock::now();

    if (argc < 4) {
        cout << "usage: " << argv[0] << " <ref.fa> <target.fa> <output.txt>" << endl;
        return 1;
    }

    cout << "using K=" << K << endl;

    compressGenome(argv[1], argv[2], argv[3]);

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "\nTotal execution time: " << elapsed.count() << " seconds\n";

    return 0;
}

//g++ -o hirgc hirgc.cpp -std=c++17
//./hirgc ref/chr22_REF.fa target/chr22_TARGET.fa chr22_output.txt
