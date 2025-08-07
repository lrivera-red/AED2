import time
import random
import os

# Tenta importar a função 'profile' para análise de memória.
# Se não for encontrado, cria uma função 'dummy' para que o script rode normalmente.
try:
    from memory_profiler import profile
except ImportError:
    def profile(func):
        return func

# -----------------------------------------------------------------------------
# ALGORITMOS DE BUSCA
# -----------------------------------------------------------------------------


def busca_forca_bruta(texto, padrao):
    n, m = len(texto), len(padrao)
    if m == 0: return []
    ocorrencias = []
    for i in range(n - m + 1):
        if texto[i:i+m] == padrao:
            ocorrencias.append(i)
    return ocorrencias

def calcular_tabela_lps(padrao):
    m = len(padrao)
    lps = [0] * m
    comprimento, i = 0, 1
    while i < m:
        if padrao[i] == padrao[comprimento]:
            comprimento += 1
            lps[i] = comprimento
            i += 1
        else:
            if comprimento != 0:
                comprimento = lps[comprimento - 1]
            else:
                lps[i] = 0
                i += 1
    return lps


def busca_kmp(texto, padrao):
    n, m = len(texto), len(padrao)
    if m == 0: return []
    lps = calcular_tabela_lps(padrao)
    ocorrencias = []
    i, j = 0, 0
    while i < n:
        if padrao[j] == texto[i]:
            i += 1
            j += 1
        if j == m:
            ocorrencias.append(i - j)
            j = lps[j - 1]
        elif i < n and padrao[j] != texto[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return ocorrencias

# -----------------------------------------------------------------------------
# FUNÇÕES UTILITÁRIAS
# -----------------------------------------------------------------------------

def ler_fasta(caminho_arquivo):
    print(f"Lendo: {caminho_arquivo}")
    try:
        with open(caminho_arquivo, 'r') as f:
            return "".join([linha.strip() for linha in f if not linha.startswith('>')]).upper()
    except FileNotFoundError:
        print(f"ERRO: Arquivo não encontrado em '{caminho_arquivo}'")
        return None

def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, '') for base in dna_sequence)[::-1]

def gerar_sequencia(tamanho, alfabeto='ACGT'):
    """Gera uma sequência aleatória."""
    return ''.join(random.choice(alfabeto) for _ in range(tamanho))

# -----------------------------------------------------------------------------
# DEFINIÇÃO DOS EXPERIMENTOS
# -----------------------------------------------------------------------------
def rodar_experimentos_sinteticos():
    """Executa um conjunto de testes com dados sintéticos."""

    print("="*60)
    print("      INICIANDO EXPERIMENTOS SINTÉTICOS EM PYTHON")
    print("="*60)

    # --- Cenário 1: Padrão fixo (m), Texto variável (n) ---
    print("\n\n--- Cenário 1: Padrão Fixo (m=100), Texto Variável (n= 10.000.000 / 20.000.000) ---")
    m_fixo = 100
    tamanhos_n = [10_000_000, 20_000_000]
    padrao_cenario1 = gerar_sequencia(m_fixo)
    for n in tamanhos_n:
        texto = gerar_sequencia(n)
        print(f"\nAnalisando n={n:,}, m={m_fixo:,}")
        
        start_time = time.time()
        busca_forca_bruta(texto, padrao_cenario1)
        print(f"  - Força Bruta: {time.time() - start_time:.4f} segundos")

        start_time = time.time()
        busca_kmp(texto, padrao_cenario1)
        print(f"  - KMP:         {time.time() - start_time:.4f} segundos")
        
    # --- Cenário 2: Texto fixo (n), Padrão variável (m) ---
    print("\n\n--- Cenário 2: Texto Fixo (n=10.000.000), Padrão Variável ---")
    n_fixo = 10_000_000
    tamanhos_m = [100, 500]
    texto_cenario2 = gerar_sequencia(n_fixo)
    for m in tamanhos_m:
        padrao = gerar_sequencia(m)
        print(f"\nAnalisando n={n_fixo:,}, m={m:,}")

        start_time = time.time()
        busca_forca_bruta(texto_cenario2, padrao)
        print(f"  - Força Bruta: {time.time() - start_time:.4f} segundos")

        start_time = time.time()
        busca_kmp(texto_cenario2, padrao)
        print(f"  - KMP:         {time.time() - start_time:.4f} segundos")

    # --- Cenário 3: Pior caso para Força Bruta ---
    print("\n\n--- Cenário 3: Pior Caso Sintético para Força Bruta (Repetições AB seguido por único C) ---")
    n, m = 20_000_000, 400
    texto_pior_caso = "AB" * (n // 2)
    padrao_pior_caso = "AB" * (m // 2 - 1) + "AC"
    print(f"\nAnalisando n={n:,}, m={m:,} (dados repetitivos)")
    
    start_time = time.time()
    busca_forca_bruta(texto_pior_caso, padrao_pior_caso)
    print(f"  - Força Bruta: {time.time() - start_time:.4f} segundos")

    start_time = time.time()
    busca_kmp(texto_pior_caso, padrao_pior_caso)
    print(f"  - KMP:         {time.time() - start_time:.4f} segundos")

def rodar_experimento_ecoli():
    print("\n\n" + "="*70)
    print("INICIANDO ANÁLISE: CASO REAL: GENOMA E. COLI VS GENE LACZ")
    print("="*70)
    
    genoma = ler_fasta("sequence_ecoli.fasta")
    padrao = ler_fasta("gene_lacz.fna")
    
    if not genoma or not padrao: return

    padrao_reverso = reverse_complement(padrao)
    print(f"Tamanho do Texto: {len(genoma):,} | Tamanho do Padrão: {len(padrao):,}")

    print("\n--- Força Bruta ---")
    start_time = time.time()
    fwd_fb = busca_forca_bruta(genoma, padrao)
    rev_fb = busca_forca_bruta(genoma, padrao_reverso)
    print(f"Tempo de Execução: {time.time() - start_time:.4f} segundos")
    print(f"Ocorrências encontradas: {len(fwd_fb)} (direta), {len(rev_fb)} (reversa)")

    print("\n--- KMP ---")
    start_time = time.time()
    fwd_kmp = busca_kmp(genoma, padrao)
    rev_kmp = busca_kmp(genoma, padrao_reverso)
    print(f"Tempo de Execução: {time.time() - start_time:.4f} segundos")
    print(f"Ocorrências encontradas: {len(fwd_kmp)} (direta), {len(rev_kmp)} (reversa)")

def rodar_experimento_humano():
    print("\n\n" + "="*70)
    print("INICIANDO ANÁLISE: CASO REAL: CROMOSSOMO 22 VS ÉXON COMT")
    print("="*70)
    
    # Baixa e descompacta o arquivo se ele não existir
    if not os.path.exists("chr22.fa"):
        print("Baixando arquivo do Cromossomo 22...")
        os.system("wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz")
        os.system("gunzip -f chr22.fa.gz")
    
    cromossomo = ler_fasta("chr22.fa")
    if not cromossomo: return
    
    padrao = "GACGCCATCACCGTGGTGACCACCAGCAACCCCAGCCTGACCGAGGACACCATCCAGGAGATGGGCCACGCCGGGGCCAAGCACGAGGGCGTGGCCGCCGACGTGGGCATCGGCCCGGAGCTGCTGGCGCCGCTGTACGACGGGCTGGGCCTGGCCAACCCCAAGGCCAAGGACATCGACACGTACGTGGAGGAGTTCTACAGCCCCCTCAAGCTC".upper()
    padrao_reverso = reverse_complement(padrao)
    print(f"Tamanho do Texto: {len(cromossomo):,} | Tamanho do Padrão: {len(padrao):,}")

    print("\n--- Força Bruta ---")
    start_time = time.time()
    fwd_fb = busca_forca_bruta(cromossomo, padrao)
    rev_fb = busca_forca_bruta(cromossomo, padrao_reverso)
    print(f"Tempo de Execução: {time.time() - start_time:.4f} segundos")
    print(f"Ocorrências encontradas: {len(fwd_fb)} (direta), {len(rev_fb)} (reversa)")

    print("\n--- KMP ---")
    start_time = time.time()
    fwd_kmp = busca_kmp(cromossomo, padrao)
    rev_kmp = busca_kmp(cromossomo, padrao_reverso)
    print(f"Tempo de Execução: {time.time() - start_time:.4f} segundos")
    print(f"Ocorrências encontradas: {len(fwd_kmp)} (direta), {len(rev_kmp)} (reversa)")

# -----------------------------------------------------------------------------
# EXECUÇÃO PRINCIPAL
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    rodar_experimento_ecoli()
    rodar_experimento_humano()
    rodar_experimentos_sinteticos()