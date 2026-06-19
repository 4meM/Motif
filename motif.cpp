#include "k-mers.h"
#include "alineamiento-multiple.h"
#include <set>
using vectorKmers =  vector<map<string, int>*>;
using posCandidato = map<string,vector<vector<int>>>;
#define ksize 7
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

vector<string> kmersCandidatos(map<string,int> analizados){
    vector<string> candidatos;
    int mayor = 0;
    for(auto it : analizados){
        if(it.second > mayor){
            mayor = it.second;
        } 
    }
    for(auto it : analizados){
        if(it.second == mayor){
            candidatos.push_back(it.first);
        }
    }
    return candidatos;
}


posCandidato findPositionCandidatos(vector<string>& candidatos, vector<string>& chains){
    posCandidato posCandidatos;
    for(int i = 0; i < candidatos.size(); i++){
        vector<vector<int>> posicionesAllSecuencias;
        string cand = candidatos[i];
        for(int j = 0; j < chains.size(); j++){
            vector<int> posUnCandidato;
            size_t pos = chains[j].find(cand, 0);
            while(pos != string::npos){
                posUnCandidato.push_back(static_cast<int>(pos));
                pos = chains[j].find(cand, pos + 1);
            }
            posicionesAllSecuencias.push_back(posUnCandidato); 
        }
        posCandidatos[cand] = posicionesAllSecuencias;
    }
    return posCandidatos;
}

void mayoresKmers(map<string,int>& analizados){
    for(auto it : analizados){
        cout << it.first << " " << it.second << "\n";
    }
    cout << "\n";
}

void printKmersCandidatos(vector<string>& candidatos){
    cout << "=========================================" << endl;
    cout << " CANDIDATOS "  << endl;
    cout << "=========================================" << endl;

    for(auto it : candidatos){
        cout<< it << endl;
    }
    cout << "\n\n\n";
}

void printPositionCandidatos(const posCandidato& misPosiciones) {
    for (const auto& parCandidato : misPosiciones) {
        cout << "=========================================" << endl;
        cout << " K-MER CANDIDATO: " << parCandidato.first << endl;
        cout << "=========================================" << endl;

        const auto& todasLasSecuencias = parCandidato.second; 
        for (size_t j = 0; j < todasLasSecuencias.size(); j++) {
            cout << "  Secuencia " << (j + 1) << " -> Posiciones: [ ";
            for (int pos : todasLasSecuencias[j]) {
                cout << pos << " ";
            }
            cout << "]" << endl;
        }
        cout << endl; 
    }
}

map<string,int> FindMaxPosComun(const posCandidato& misPosiciones) {
    map<string,int> maxPosComun;
    for(const auto pos: misPosiciones){
        auto v=pos.second;
        map<int,int> frecuencia;
        for(int i=0;i<v.size();i++){
            for(int j=0;j<v[i].size();j++){
                frecuencia[v[i][j]]++;
            }
        }
        int maxFrec=0;
        for(auto i:frecuencia){
            if(maxFrec<i.second){
                maxPosComun[pos.first]=i.first;
                maxFrec=i.second;
            }
        }
    }
    return maxPosComun;
}
void mostrarPosFrec(const map<string, int>& maxPosComun) {
    cout << "===== POSICION MAS FRECUENTE POR ELEMENTO =====\n\n";

    for (const auto& elemento : maxPosComun) {
        cout << "Clave: " << elemento.first << '\n';
        cout << "Posicion mas comun: " << elemento.second << '\n';
        cout << "----------------------------------------\n";
    }
}

map<string,vector<string>> ExtraerRegCommon(vector<string> chains, int k, map<string,int> startSecuence){
    map<string,vector<string>> RegCommonExt;
    for(auto pos:startSecuence){
        size_t inicio=pos.second;
        vector<string> auxsec;
        for(string sec: chains){
            string subcadena= sec.substr(inicio,k);
            auxsec.push_back(subcadena);
        }
        RegCommonExt[pos.first]=auxsec;
    }
    return RegCommonExt;
}
void mostrarRegComExtr(const map<string, vector<string>>& RegCommonExt) {
    cout << "===== REGIONES COMUNES EXTRAIDAS =====\n\n";

    for (const auto& elemento : RegCommonExt) {
        cout << "Clave: " << elemento.first << '\n';
        cout << "Subcadenas extraidas:\n";

        int i = 1;
        for (const string& subcadena : elemento.second) {
            cout << "   " << i++ << ". " << subcadena << '\n';
        }

        cout << "----------------------------------------\n";
    }
}

map<string, int> alineamientoMultple(map<string,vector<string>>& RegCommonExt){
    map<string, int> resultado; // Mapa para retornar al final

    for(auto& it : RegCommonExt){ 
        
        TreeGuide(it.second);
        AlineamientoMultiple();
        
        for(const auto& entry : clustersMap.begin()->first->chains){
            cout << "Cadena alineada: " << *entry << "\n";
        } 

        clustersMap.clear(); 
    }

    return resultado; 
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
    allKmers(vk,chains,ksize);
    map<string,int> analizados = kmersFrecuente(vk);
    mayoresKmers(analizados);
    vector<string> candidatos = kmersCandidatos(analizados);
    printKmersCandidatos(candidatos);
    posCandidato posCandidatos = findPositionCandidatos(candidatos, chains);
    printPositionCandidatos(posCandidatos);
    map<string,int> MaxposKmer=FindMaxPosComun(posCandidatos);
    mostrarPosFrec(MaxposKmer);
    map<string,vector<string>> RegCommExtr = ExtraerRegCommon(chains,ksize,MaxposKmer);
    mostrarRegComExtr(RegCommExtr);
    alineamientoMultple(RegCommExtr);
}