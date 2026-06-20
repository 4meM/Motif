#include "k-mers.h"
#include "alineamiento-multiple.h"
#include <set>
using vectorKmers =  vector<map<string, int>*>;
using posCandidato = map<string,vector<vector<int>>>;
#define ksize 9
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

vector<map<char,int>> construirMatrizFrecuencias(const vector<string>& regiones) {
    if (regiones.empty()) return {};
    int largo = regiones[0].size();
    vector<map<char,int>> matriz(largo);
    for (const string& seq : regiones) {
        for (int i = 0; i < largo; ++i) {
            matriz[i][seq[i]]++;
        }
    }
    return matriz;
}

string obtenerConsenso(const vector<map<char,int>>& matriz) {
    string consenso;
    for (const auto& posMap : matriz) {
        char baseMax = 'N';
        int maxCount = 0;
        for (const auto& par : posMap) {
            if (par.second > maxCount) {
                maxCount = par.second;
                baseMax = par.first;
            }
        }
        consenso += baseMax;
    }
    return consenso;
}

void imprimirMatrizFrecuencias(const vector<map<char,int>>& matriz) {
    if (matriz.empty()) {
        cout << "La matriz está vacía." << endl;
        return;
    }

    cout << "\n===== MATRIZ DE FRECUENCIAS =====\n";
    cout << "Posicion |  A  |  C  |  G  |  T  | Consenso\n";
    cout << "---------|-----|-----|-----|-----|---------\n";

    for (size_t i = 0; i < matriz.size(); ++i) {
        int countA = 0, countC = 0, countG = 0, countT = 0;
        for (const auto& par : matriz[i]) {
            if (par.first == 'A') countA = par.second;
            else if (par.first == 'C') countC = par.second;
            else if (par.first == 'G') countG = par.second;
            else if (par.first == 'T') countT = par.second;
        }

        char consenso = 'N';
        int maxCount = 0;
        if (countA > maxCount) { maxCount = countA; consenso = 'A'; }
        if (countC > maxCount) { maxCount = countC; consenso = 'C'; }
        if (countG > maxCount) { maxCount = countG; consenso = 'G'; }
        if (countT > maxCount) { maxCount = countT; consenso = 'T'; }

        cout << "    " << i+1 << "    |";
        cout << "  " << countA << "  |";
        cout << "  " << countC << "  |";
        cout << "  " << countG << "  |";
        cout << "  " << countT << "  |";
        cout << "    " << consenso << "\n";
    }
    cout << "-----------------------------------------\n";
    cout << "Total de secuencias: " << matriz[0].size() << "\n\n";
}

vector<double> calcularPorcentajesConservacion(const vector<map<char,int>>& matriz, int totalSecuencias) {
    vector<double> porcentajes;
    for (const auto& posMap : matriz) {
        int maxCount = 0;
        for (const auto& par : posMap) {
            if (par.second > maxCount) maxCount = par.second;
        }
        double pct = (totalSecuencias > 0) ? (maxCount * 100.0 / totalSecuencias) : 0.0;
        porcentajes.push_back(pct);
    }
    return porcentajes;
}

void identificarPosicionesConservadasVariables(const vector<map<char,int>>& matriz, int totalSecuencias,vector<int>& posConservadas,vector<int>& posVariables) {
    posConservadas.clear();
    posVariables.clear();
    for (size_t i = 0; i < matriz.size(); ++i) {
        int maxCount = 0;
        for (const auto& par : matriz[i]) {
            if (par.second > maxCount) maxCount = par.second;
        }
        if (maxCount == totalSecuencias) {
            posConservadas.push_back(i);
        } else {
            posVariables.push_back(i);
        }
    }
}

void describirBasesVariables(const vector<map<char,int>>& matriz, const vector<int>& posVariables) {
    cout << "===== BASES OBSERVADAS EN POSICIONES VARIABLES =====\n";
    for (int pos : posVariables) {
        cout << "Posicion " << pos+1 << ": ";
        for (const auto& par : matriz[pos]) {
            cout << par.first << "(" << par.second << ") ";
        }
        cout << endl;
    }
    cout << endl;
}

double porcentajeGlobalConservacion(const vector<double>& porcentajes) {
    if (porcentajes.empty()) return 0.0;
    double suma = 0.0;
    for (double p : porcentajes) suma += p;
    return suma / porcentajes.size();
}

void reportarMotif(const string& consenso, int longitud, 
                   const map<string,int>& posInicioPorSecuencia,
                   int totalSecuenciasConMotif,
                   const vector<double>& porcentajes) {
    cout << "\n===== REPORTE DEL MOTIF =====\n";
    cout << "Motif (secuencia consenso): " << consenso << endl;
    cout << "Longitud del motif: " << longitud << endl;
    cout << "Numero de secuencias que contienen el motif: " << totalSecuenciasConMotif << endl;
    cout << "Posiciones de inicio en cada secuencia:\n";
    for (const auto& par : posInicioPorSecuencia) {
        cout << "  " << par.first << ": " << par.second << endl;
    }
    double global = porcentajeGlobalConservacion(porcentajes);
    cout << "Porcentaje global de conservacion: " << global << "%" << endl;
    cout << "================================\n";
}

int main(){
    /*vector<string> chains = {
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
    };*/
    vector<string> chains = {
        "ATCGGCTAACGTAGCTAGCTTGACCGTACGATCGATCGGATCGTAGCTAGCATCGATCGTACGATCGATGCTAGCTAGCATCGATCGATACGATCGTAGCTAGCTACGTAGCTAGCTACGTAGCTTACGATGACGGTACCGATCGATCGTAGCTAACGTA",
        "GCTAGCTAGCATCGATCGTAGCTAGCTAGCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATACGTAGCTACGTACGATGACATCGTAGCTAGCTAACGTAGCTAGCTAGCGATCGTAGCTAGCTA",
        "CGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTACGTACGATGATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTA",
        "AACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTACGTATACGATGACGCTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGT",
        "TAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAACGTATACGATGTCCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACG",
        "GATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTATACGATGACGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCT",
        "CTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTATACGATGCCGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTA",
        "CGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAACGTAGCTAGCTAGCATCGATCGTAGCTATACGATGACCGTAGCTAACGTAGCTAGCATCGATCGTAGCTAAC",
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
    
    string primerKmer = RegCommExtr.begin()->first;
    vector<string> regiones = RegCommExtr[primerKmer];

    vector<map<char,int>> matriz = construirMatrizFrecuencias(regiones);
    imprimirMatrizFrecuencias(matriz);
    string consenso = obtenerConsenso(matriz);
    cout << "\nSecuencia consenso: " << consenso << endl;

    int totalSecuencias = regiones.size();
    vector<double> porcentajes = calcularPorcentajesConservacion(matriz, totalSecuencias);
    cout << "\nPorcentajes de conservacion por posicion:\n";
    for (size_t i = 0; i < porcentajes.size(); ++i) {
        cout << "  Pos " << i+1 << ": " << porcentajes[i] << "%" << endl;
    }

    vector<int> posConservadas, posVariables;
    identificarPosicionesConservadasVariables(matriz, totalSecuencias, posConservadas, posVariables);
    cout << "\nPosiciones totalmente conservadas (100%): ";
    for (int p : posConservadas) cout << p+1 << " ";
    cout << endl;
    cout << "Posiciones variables: ";
    for (int p : posVariables) cout << p+1 << " ";
    cout << endl;

    describirBasesVariables(matriz, posVariables);

    double global = porcentajeGlobalConservacion(porcentajes);
    cout << "\nPorcentaje global de conservacion: " << global << "%" << endl;

    map<string,int> posInicioPorSecuencia;
    for (size_t i = 0; i < chains.size(); ++i) {
        string secName = "S" + to_string(i+1);
        posInicioPorSecuencia[secName] = MaxposKmer[primerKmer];
    }
    reportarMotif(consenso, consenso.size(), posInicioPorSecuencia, totalSecuencias, porcentajes);
}