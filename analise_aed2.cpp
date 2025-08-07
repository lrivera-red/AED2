#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <random>

// -----------------------------------------------------------------------------
// ALGORITMOS DE BUSCA
// -----------------------------------------------------------------------------

std::vector<int> busca_forca_bruta(const std::string& texto, const std::string& padrao) {
    int n = texto.length();
    int m = padrao.length();
    std::vector<int> ocorrencias;
    if (m == 0) return ocorrencias;

    for (int i = 0; i <= n - m; ++i) {
        if (texto.substr(i, m) == padrao) {
            ocorrencias.push_back(i);
        }
    }
    return ocorrencias;
}

std::vector<int> calcular_tabela_lps(const std::string& padrao) {
    int m = padrao.length();
    std::vector<int> lps(m, 0);
    int comprimento = 0;
    int i = 1;
    while (i < m) {
        if (padrao[i] == padrao[comprimento]) {
            comprimento++;
            lps[i] = comprimento;
            i++;
        } else {
            if (comprimento != 0) {
                comprimento = lps[comprimento - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
    return lps;
}

std::vector<int> busca_kmp(const std::string& texto, const std::string& padrao) {
    int n = texto.length();
    int m = padrao.length();
    std::vector<int> ocorrencias;
    if (m == 0) return ocorrencias;

    std::vector<int> lps = calcular_tabela_lps(padrao);
    int i = 0, j = 0;
    while (i < n) {
        if (padrao[j] == texto[i]) {
            i++;
            j++;
        }
        if (j == m) {
            ocorrencias.push_back(i - j);
            j = lps[j - 1];
        } else if (i < n && padrao[j] != texto[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    return ocorrencias;
}

// -----------------------------------------------------------------------------
// FUNÇÕES UTILITÁRIAS
// -----------------------------------------------------------------------------

std::string ler_fasta(const std::string& caminho_arquivo) {
    std::cout << "Lendo: " << caminho_arquivo << std::endl;
    std::ifstream arquivo(caminho_arquivo);
    if (!arquivo) {
        std::cerr << "ERRO: Arquivo não encontrado em '" << caminho_arquivo << "'" << std::endl;
        return "";
    }
    std::stringstream buffer;
    std::string linha;
    while (std::getline(arquivo, linha)) {
        if (linha.empty() || linha[0] != '>') {
            buffer << linha;
        }
    }
    std::string sequencia = buffer.str();
    std::transform(sequencia.begin(), sequencia.end(), sequencia.begin(), ::toupper);
    return sequencia;
}

std::string gerar_sequencia(int tamanho, const std::string& alfabeto) {
    std::string resultado;
    resultado.reserve(tamanho);
    std::random_device rd;
    std::mt19937 gerador(rd());
    std::uniform_int_distribution<> distribuicao(0, alfabeto.length() - 1);
    for (int i = 0; i < tamanho; ++i) {
        resultado += alfabeto[distribuicao(gerador)];
    }
    return resultado;
}

std::string reverse_complement(const std::string& dna) {
    std::string reverso;
    reverso.reserve(dna.length());
    for (int i = dna.length() - 1; i >= 0; --i) {
        switch (dna[i]) {
            case 'A': reverso += 'T'; break;
            case 'T': reverso += 'A'; break;
            case 'C': reverso += 'G'; break;
            case 'G': reverso += 'C'; break;
            default:  reverso += 'N'; break;
        }
    }
    return reverso;
}

// -----------------------------------------------------------------------------
// DEFINIÇÃO DOS EXPERIMENTOS
// -----------------------------------------------------------------------------

void rodar_experimentos_sinteticos() {
    std::cout << std::fixed << std::setprecision(4);
    
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "      INICIANDO EXPERIMENTOS SINTÉTICOS EM C++" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    // --- Cenário 1: Padrão fixo (m), Texto variável (n) ---
    std::cout << "\n\n--- Cenário 1: Padrão Fixo (m=100), Texto Variável ---" << std::endl;
    int m_fixo = 100;
    std::vector<long long> tamanhos_n = {10000000, 20000000};
    std::string padrao_cenario1 = gerar_sequencia(m_fixo, "ACGT");
    for (long long n : tamanhos_n) {
        std::string texto = gerar_sequencia(n, "ACGT");
        std::cout << "\nAnalisando n=" << n << ", m=" << m_fixo << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        busca_forca_bruta(texto, padrao_cenario1);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> t = end - start;
        std::cout << "  - Força Bruta: " << t.count() << " segundos" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        busca_kmp(texto, padrao_cenario1);
        end = std::chrono::high_resolution_clock::now();
        t = end - start;
        std::cout << "  - KMP:         " << t.count() << " segundos" << std::endl;
    }

    // --- Cenário 2: Texto fixo (n), Padrão variável (m) ---
    std::cout << "\n\n--- Cenário 2: Texto Fixo (n=10,000,000), Padrão Variável ---" << std::endl;
    long long n_fixo = 10000000;
    std::vector<int> tamanhos_m = {100, 500};
    std::string texto_cenario2 = gerar_sequencia(n_fixo, "ACGT");
    for (int m : tamanhos_m) {
        std::string padrao = gerar_sequencia(m, "ACGT");
        std::cout << "\nAnalisando n=" << n_fixo << ", m=" << m << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        busca_forca_bruta(texto_cenario2, padrao);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> t = end - start;
        std::cout << "  - Força Bruta: " << t.count() << " segundos" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        busca_kmp(texto_cenario2, padrao);
        end = std::chrono::high_resolution_clock::now();
        t = end - start;
        std::cout << "  - KMP:         " << t.count() << " segundos" << std::endl;
    }

    // --- Cenário 3: Pior caso para Força Bruta ---
    std::cout << "\n\n--- Cenário 3: Pior Caso Sintético para Força Bruta ---" << std::endl;
    long long n_pior = 20000000;
    int m_pior = 400;
    
    std::string texto_pior_caso;
    texto_pior_caso.reserve(n_pior);
    for (long long i = 0; i < n_pior / 2; ++i) texto_pior_caso += "AB";
    
    std::string padrao_pior_caso;
    padrao_pior_caso.reserve(m_pior);
    for (int i = 0; i < m_pior / 2 - 1; ++i) padrao_pior_caso += "AB";
    padrao_pior_caso += "AC";

    std::cout << "\nAnalisando n=" << n_pior << ", m=" << m_pior << " (dados repetitivos)" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    busca_forca_bruta(texto_pior_caso, padrao_pior_caso);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t = end - start;
    std::cout << "  - Força Bruta: " << t.count() << " segundos" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    busca_kmp(texto_pior_caso, padrao_pior_caso);
    end = std::chrono::high_resolution_clock::now();
    t = end - start;
    std::cout << "  - KMP:         " << t.count() << " segundos" << std::endl;
}

void rodar_experimento_ecoli() {
    std::cout << "\n\n" << std::string(70, '=') << std::endl;
    std::cout << "INICIANDO ANÁLISE: CASO REAL: GENOMA E. COLI VS GENE LACZ" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::string genoma = ler_fasta("sequence_ecoli.fasta");
    std::string padrao = ler_fasta("gene_lacz.fna");
    
    if (genoma.empty() || padrao.empty()) return;

    std::string padrao_reverso = reverse_complement(padrao);

    std::cout << "Tamanho do Texto: " << genoma.length() << " | Tamanho do Padrão: " << padrao.length() << std::endl;

    std::cout << "\n--- Força Bruta ---" << std::endl;
    auto start_fb = std::chrono::high_resolution_clock::now();
    auto fwd_fb = busca_forca_bruta(genoma, padrao);
    auto rev_fb = busca_forca_bruta(genoma, padrao_reverso);
    auto end_fb = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_fb = end_fb - start_fb;
    std::cout << "Tempo de Execução: " << std::fixed << std::setprecision(4) << t_fb.count() << " segundos" << std::endl;
    std::cout << "Ocorrências encontradas: " << fwd_fb.size() << " (direta), " << rev_fb.size() << " (reversa)" << std::endl;

    std::cout << "\n--- KMP ---" << std::endl;
    auto start_kmp = std::chrono::high_resolution_clock::now();
    auto fwd_kmp = busca_kmp(genoma, padrao);
    auto rev_kmp = busca_kmp(genoma, padrao_reverso);
    auto end_kmp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_kmp = end_kmp - start_kmp;
    std::cout << "Tempo de Execução: " << std::fixed << std::setprecision(4) << t_kmp.count() << " segundos" << std::endl;
    std::cout << "Ocorrências encontradas: " << fwd_kmp.size() << " (direta), " << rev_kmp.size() << " (reversa)" << std::endl;
}

void rodar_experimento_humano() {
    std::cout << "\n\n" << std::string(70, '=') << std::endl;
    std::cout << "INICIANDO ANÁLISE: CASO REAL: CROMOSSOMO 22 VS ÉXON COMT" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Baixa e descompacta o arquivo do Cromossomo 22
    system("wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz");
    system("gunzip -f chr22.fa.gz");
    
    std::string cromossomo = ler_fasta("chr22.fa");
    if (cromossomo.empty()) return;
    
    std::string padrao = "GACGCCATCACCGTGGTGACCACCAGCAACCCCAGCCTGACCGAGGACACCATCCAGGAGATGGGCCACGCCGGGGCCAAGCACGAGGGCGTGGCCGCCGACGTGGGCATCGGCCCGGAGCTGCTGGCGCCGCTGTACGACGGGCTGGGCCTGGCCAACCCCAAGGCCAAGGACATCGACACGTACGTGGAGGAGTTCTACAGCCCCCTCAAGCTC";
    std::transform(padrao.begin(), padrao.end(), padrao.begin(), ::toupper);
    std::string padrao_reverso = reverse_complement(padrao);

    std::cout << "Tamanho do Texto: " << cromossomo.length() << " | Tamanho do Padrão: " << padrao.length() << std::endl;

    std::cout << "\n--- Força Bruta ---" << std::endl;
    auto start_fb = std::chrono::high_resolution_clock::now();
    auto fwd_fb = busca_forca_bruta(cromossomo, padrao);
    auto rev_fb = busca_forca_bruta(cromossomo, padrao_reverso);
    auto end_fb = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_fb = end_fb - start_fb;
    std::cout << "Tempo de Execução: " << std::fixed << std::setprecision(4) << t_fb.count() << " segundos" << std::endl;
    std::cout << "Ocorrências encontradas: " << fwd_fb.size() << " (direta), " << rev_fb.size() << " (reversa)" << std::endl;

    std::cout << "\n--- KMP ---" << std::endl;
    auto start_kmp = std::chrono::high_resolution_clock::now();
    auto fwd_kmp = busca_kmp(cromossomo, padrao);
    auto rev_kmp = busca_kmp(cromossomo, padrao_reverso);
    auto end_kmp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_kmp = end_kmp - start_kmp;
    std::cout << "Tempo de Execução: " << std::fixed << std::setprecision(4) << t_kmp.count() << " segundos" << std::endl;
    std::cout << "Ocorrências encontradas: " << fwd_kmp.size() << " (direta), " << rev_kmp.size() << " (reversa)" << std::endl;
}

// -----------------------------------------------------------------------------
// EXECUÇÃO PRINCIPAL
// -----------------------------------------------------------------------------

int main() {
    rodar_experimento_ecoli();
    rodar_experimento_humano();
    rodar_experimentos_sinteticos();

    return 0;
}