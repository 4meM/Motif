#ifndef K_MERS_H
#define K_MERS_H

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <set>
#include <chrono>

using namespace std;

void topN(int N, const map<string, int>& freqs) {
    multimap<int, string> maxMin;
    for (const auto& pair : freqs) {
        maxMin.insert({pair.second, pair.first});
    }
    cout << "\n--- Top " << N << " K-mers ---\n";
    int contador = 0;
    for (auto it = maxMin.rbegin(); it != maxMin.rend() && contador < N; ++it, ++contador) {
        cout << it->second << ": " << it->first << "\n";
    }
    cout << "-----------------------\n";
}

string readFastaFile(const string& archivo, string& dna_out) {
    dna_out.clear();
    ifstream fastaFile(archivo);
    if (!fastaFile.is_open()) {
        cerr << "Error al abrir el archivo: " << archivo << endl;
        return "";
    }
    string line, id_file = "";
    bool id_found = false;
    while (getline(fastaFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!id_found) {
                size_t first_space = line.find(' ');
                id_file = line.substr(1, first_space != string::npos ? first_space - 1 : string::npos);
                id_found = true;
            }
        } else {
            dna_out += line;
        }
    }
    return id_file;
}

void conteo(const string& dna, int k, map<string, int>& freqs_out) {
    freqs_out.clear();
    if (dna.length() < k) return;
    for (size_t i = 0; i <= dna.length() - k; ++i) {
        freqs_out[dna.substr(i, k)]++;
    }
}

void show_info(const string& archivo, int k) {
    cout << "\n================================================\n";
    cout << "INFORME INDIVIDUAL PARA: " << archivo << " (k=" << k << ")\n";
    cout << "================================================\n";
    
    string dna;
    auto inicio_lectura = chrono::high_resolution_clock::now();
    string id = readFastaFile(archivo, dna);
    
    map<string, int> freqs;
    conteo(dna, k, freqs);
    auto fin_conteo = chrono::high_resolution_clock::now();
    auto duracion = chrono::duration_cast<chrono::milliseconds>(fin_conteo - inicio_lectura);

    cout << "ID de la secuencia: " << id << endl;
    cout << "Numero de K-mers distintos: " << freqs.size() << "\n";
    cout << "DURACION (Lectura y Conteo): " << duracion.count() << " ms\n";
    
    topN(20, freqs);
}


#endif
/*
int main() {
    show_info("compare-fasta-files/sequence.fasta", 21);
    //show_info("compare-fasta-files/sequence2.fasta", 21);
    int k = 21;
    vector<string> archivos_a_comparar = {
        "compare-fasta-files/sequence.fasta",
        "compare-fasta-files/sequence2.fasta",
    };

    vector<map<string, int>> perfiles_de_frecuencia;
    vector<string> ids_de_secuencias;

    cout << "Generando perfiles para comparacion (k=" << k << ")..." << endl;
    for (const string& archivo : archivos_a_comparar) {
        string dna;
        string id = readFastaFile(archivo, dna);
        map<string, int> freqs;
        conteo(dna, k, freqs);
        perfiles_de_frecuencia.push_back(freqs);
        ids_de_secuencias.push_back(id);
        cout << "  - Perfil generado para '" << id << "'" << endl;
    }

    realizarAnalisisComparativo(perfiles_de_frecuencia, ids_de_secuencias);
    return 0;
}
*/