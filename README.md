# Bioinformatika
Genome compression and decompression

Testni primjeri i rezultati se mogu preuzeti sa linka: 
https://drive.google.com/file/d/1Nnk5JMTWkxEz3OhWI5vRXBI8DRa5vbxa/view?usp=sharing

Pokretanje:
1. Kompresija
g++ -o hirgc hirgc.cpp -std=c++17
./hirgc ref/reference_file.fa target/target_file.fa output.txt
2. Dekompresija
g++ -o hirgdc hirgdc.cpp -std=c++17
./hirgdc.exe output.txt.7z recovered.fa
