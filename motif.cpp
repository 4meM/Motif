#include "k-mers.h"
#include "alineamiento-multiple.h"
#include <set>
using vectorKmers =  vector<map<string, int>*>;


void allKmers(vectorKmers &perfiles, vector<string>& chains, int k){
    for(int i = 0; i < chains.size(); i++){
        map<string, int>* perfil = new map<string, int>();
        conteo(chains[i], k, *perfil);
        perfiles.push_back(perfil);
    }
}

map<string, int> kmersFrecuente(const vectorKmers &perfiles) {
    map<string, int> analizados; 

    for (size_t i = 0; i < perfiles.size(); i++) {
        map<string, int>* mapaActual = perfiles[i];

        for (const auto& pairSI : *mapaActual) {
            string kmer = pairSI.first;

            if (analizados.find(kmer) != analizados.end()) {
                continue; 
            }

            analizados[kmer] = 1;

            for (size_t j = 0; j < perfiles.size(); j++) {
                if (i == j) continue; 
                map<string, int>* otroMapa = perfiles[j];

                if (otroMapa->find(kmer) != otroMapa->end()) analizados[kmer]++;
            }
        }
    }

    return analizados; 
}

void mayoresKmers(map<string,int>& analizados){
    for(auto it : analizados){
        cout << it.first << " " << it.second << "\n";
    }
    cout << "\n";
}




int main(){
    vector<string> chains = {
        "ATCGTACGATGACCTGATCG",
        "GGTATACGATGACGTTACCA",
        "TTTCTACGATGACCATAGGT",
        "AACGTACGATGACGGGTTAA",
        "CGGATACGATGACTTCCGTA",
        "TACCTACGATGACAGGTACA",
        "GACTTACGATGACCGATAGC",
        "TCGATACGATGACTGGCAAT",
        "AGGCTACGATGACATTCGGA",
        "CCTATACGATGACGGAATTC",
    };
    vectorKmers vk;
    allKmers(vk,chains,7);
    map<string,int> analizados = kmersFrecuente(vk);
    mayoresKmers(analizados);
}