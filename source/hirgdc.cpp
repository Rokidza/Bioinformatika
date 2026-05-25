//g++ -o hirgc hirgc.cpp -std=c++17
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <chrono>

using namespace std;

// result = "AGCGATCGATCATGC.."
string readGenomeFile(string filename) {
    ifstream f(filename);
    if (!f.is_open()) {
        cout << "ERROR: could not open file: " << filename << endl;
        exit(1);
    }

    string result = "";
    string line = "";
    while (getline(f, line)) {
        if (line.length() == 0) continue;
        if (line[0] == '>') continue;
        for (int i = 0; i < (int)line.length(); i++) {
            line[i] = toupper(line[i]);
        }
        result += line;
    }
    f.close();
    return result;
}

void decompressGenome(string inputFile, string outputFile) {
    string txtFile = inputFile;

    // extractaj .7z file u txt file
    if (inputFile.length() > 3 && inputFile.substr(inputFile.length() - 3) == ".7z") {
        cout << "Extracting from 7-zip file..." << endl;

        string cmd = "7z e -y " + inputFile;
        int result = system(cmd.c_str());

        if (result != 0) {
            cout << "7-Zip extraction failed." << endl; 
        }
        txtFile = inputFile.substr(0, inputFile.length() - 3); //file.txt bez .7z
    }

    ifstream inFile(txtFile);
    if (!inFile.is_open()) {
        cout << "ERROR: cant open text file " << txtFile << endl;
        exit(1);
    }

    string magic;
    inFile >> magic;
    if (magic != "HIRG") { 
        cout << "ERROR: Not a valid text HIRG file, or corrupted." << endl;
        exit(1);
    }
    
    // nepotreban ostatak prvog reda
    string dummy;
    getline(inFile, dummy);
    
    // 2. redak -> metadata
    string targetMetadata;
    getline(inFile, targetMetadata);

    // isčita sve Header info koje smo zapisali
    string refFile;
    inFile >> refFile;

    long long targetLen;
    inFile >> targetLen;
    
    int kValue;
    inFile >> kValue;
    
    size_t numBlocks;
    inFile >> numBlocks;
    
    struct Block {
        size_t pos;
        size_t len;
    };
    vector<Block> nBlocks(numBlocks);
    //spremi pos i len za svaki N blok.
    for (size_t i = 0; i < numBlocks; i++) {
        inFile >> nBlocks[i].pos >> nBlocks[i].len;
    }

    size_t numLowerBlocks = 0;
    inFile >> numLowerBlocks;
    vector<Block> lowerBlocks(numLowerBlocks);
    // za restoreat upper i lowercase slova
    for (size_t i = 0; i < numLowerBlocks; i++) {
        inFile >> lowerBlocks[i].pos >> lowerBlocks[i].len;
    }

    cout << "Reference file : " << refFile << endl;
    cout << "Target length  : " << targetLen << endl;
    cout << "K value used   : " << kValue << endl;
    cout << "N blocks found : " << nBlocks.size() << endl;

    cout << "Loading reference genome..." << endl;
    string ref = readGenomeFile(refFile);
    
    cout << "Decompressing..." << endl;
    string recovered;
    recovered.reserve(targetLen);

    long long pos = 0;
    while (pos < targetLen) {
        char tag; // u tag spremamo M za match, L za literal
        if (!(inFile >> tag)) break;
        
        if (tag == 'M') {
            uint32_t refPos;
            size_t matchLen;
            inFile >> refPos >> matchLen;
            
            recovered.append(ref, refPos, matchLen); // appenda tocno taj dio iz ref
            pos += matchLen;
        } else if (tag == 'L') {
            char c;
            inFile >> c;
            recovered.push_back(c); // doda jedan literal iz {A,C,G,T}
            pos++;
        }
    }
    
    inFile.close();
    
    // recovera N znakove
    for (const auto& block : nBlocks) {
        for (size_t i = 0; i < block.len; i++) {
            if (block.pos + i < recovered.length()) {
                recovered[block.pos + i] = 'N';
            }
        }
    }

    // recovera lower case znakove
    for (const auto& block : lowerBlocks) {
        for (size_t i = 0; i < block.len; i++) {
            if (block.pos + i < recovered.length()) {
                recovered[block.pos + i] = tolower(recovered[block.pos + i]);
            }
        }
    }
    
    cout << "Writing recovered genome to " << outputFile << "..." << endl;
    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cout << "ERROR: cant open output file " << outputFile << endl;
        exit(1);
    }
    

    // zapiši u output file
    if (targetMetadata.empty() || targetMetadata[0] != '>') {
        outFile << ">" << targetMetadata << "\n";
    } else {
        outFile << targetMetadata << "\n";
    }
    for (size_t i = 0; i < recovered.length(); i += 50) {
        outFile << recovered.substr(i, 50) << "\n";     // u svaki redak zapise 50 slova
    }
    outFile.close();
    
    cout << "Decompression done." << endl;
}

int main(int argc, char* argv[]) {
    auto start_time = chrono::high_resolution_clock::now();

    if (argc < 3) {
        cout << "usage example: " << argv[0] << " out.txt.7z recovered.fa" << endl;
        return 1;
    }

    decompressGenome(argv[1], argv[2]);

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "\nTotal execution time: " << elapsed.count() << " seconds\n";

    return 0;
}

//g++ -o hirgc hirgc.cpp -std=c++17
