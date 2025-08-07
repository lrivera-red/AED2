# Análise Comparativa de Algoritmos de Busca: Força Bruta vs. KMP

![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![C++](https://img.shields.io/badge/C%2B%2B-00599C?style=for-the-badge&logo=c%2B%2B&logoColor=white)

## Descrição do Projeto

Este repositório contém os códigos-fonte para o trabalho da disciplina de Algoritmos e Estruturas de Dados II, que realiza uma análise comparativa teórica e prática entre os algoritmos de busca de padrões **Força Bruta** e **Knuth-Morris-Pratt (KMP)**.

O estudo de caso foca em uma aplicação realística de bioinformática: a busca por sequências genéticas (padrões) em genomas (textos). O objetivo é determinar qual algoritmo é mais efetivo, considerando não apenas a complexidade assintótica, mas também o tempo de execução prático e o impacto do ambiente de implementação.

A análise é conduzida em duas linguagens distintas:
1.  **Python 3:** Uma linguagem interpretada de alto nível, popular em bioinformática.
2.  **C++17:** Uma linguagem compilada de baixo nível, utilizada para validar a performance teórica dos algoritmos.

A principal descoberta do trabalho é a notável divergência entre a performance teórica e a prática, onde o Força Bruta se mostra surpreendentemente mais rápido em Python, enquanto o KMP demonstra sua superioridade massiva em C++.

## Funcionalidades

-   Implementação dos algoritmos Força Bruta e Knuth-Morris-Pratt (KMP).
-   Análise de performance (tempo de execução) em Python e C++.
-   Testes com dados sintéticos para validar casos de borda (pior caso).
-   Testes com dados genômicos reais (genoma da *E. coli* e Cromossomo 22 Humano).
-   Script Python inclui análise de consumo de memória.

## Pré-requisitos

Para executar os experimentos, o seguinte ambiente é necessário (testado em **Ubuntu 22.04**):

#### Ambiente Geral
-   `git` para clonar o repositório.
-   `wget` e `gunzip` para baixar e descompactar os dados do genoma humano (geralmente já vêm instalados).

#### Para o código Python
-   Python 3.8 ou superior.
-   `pip` (gerenciador de pacotes do Python).
-   A biblioteca `memory-profiler` para análise de espaço.

```bash
# Instalar a biblioteca de perfil de memória
pip install memory-profiler
```

#### Para o código C++
-   Compilador C++ com suporte a C++17 (como o `g++`).

```bash
# Instalar o compilador g++
sudo apt update
sudo apt install build-essential
```

## Como Executar

### 1. Clone o Repositório

```bash
git clone <URL_DO_SEU_REPOSITORIO>
cd <NOME_DA_PASTA_DO_REPOSITORIO>
```

### 2. Prepare os Arquivos de Dados

Para o teste com o genoma da *E. coli*, posicione os arquivos `sequence_ecoli.fasta` e `gene_lacz.fna` dentro da pasta do projeto. O arquivo do Cromossomo 22 (`chr22.fa`) será baixado automaticamente pelos scripts.

### 3. Executando a Análise em C++

Primeiro, compile o código-fonte. A flag `-O3` é recomendada para otimizar o código para máxima velocidade.

```bash
# Compilar o código
g++ -std=c++17 -O3 -o analise_aed2 analise_aed2.cpp

# Executar a análise
./analise_aed2
```
O programa executará os três cenários de teste (E. coli, Cromossomo 22, e Sintético) e imprimirá os tempos de execução no terminal.

### 4. Executando a Análise em Python

O script Python pode ser executado de duas formas:

**A) Para Análise de Tempo de Execução:**
```bash
python3 analise_aed2.py
```
Este comando executará todos os cenários de teste e imprimirá apenas os resultados de tempo de execução.

**B) Para Análise de Tempo e Espaço (Memória):**
```bash
python3 -m memory_profiler analise_aed2.py
```
Este comando executará o script através do `memory-profiler`. Ao final da execução normal, ele imprimirá um relatório detalhado do uso de memória das funções `busca_forca_bruta` e `busca_kmp`.

## Estrutura do Repositório

```
/
|-- analise_aed2.py         # Script principal para análise em Python
|-- analise_aed2.cpp        # Código-fonte principal para análise em C++
|-- sequence_ecoli.fasta    # [Necessário] Dados do genoma da E. coli (texto)
|-- gene_lacz.fna           # [Necessário] Dados do gene lacZ (padrão)
|-- chr22.fa                # (Baixado automaticamente pelos scripts)
|-- README.md               # Este arquivo
```

## Autor
-   Luiz Gabriel Silva Rivera
