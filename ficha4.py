# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 18:29:17 2020

@author: Dias
"""


#%%
import re
DNA = 'ACTG'
RNA = 'ACUG'
PROT = 'ACDEFGHIKLMNPQRSTVWY_'


"""Função que valida a sequencia colocada tendo como argumentos uma string com a sequencia
e o segundo argumento o tipo de sequencia que se quer procurar validar (DNA, RNA ou PROT).
Returna o booleano True caso se verificar que seja valido e False caso o contrário.
"""
def valida(seq,alf):
    seq = seq.upper()
    res = True
    i = 0
    while i < len(seq) and res:
        if seq[i] not in alf:
            res = False
        else: i += 1
    return res


#%%
"""
Função que faz o return do tipo de sequencia que se está a analisar.
"""
def alfab(seq):
    seq = seq.upper()
    if valida(seq,DNA):
        return DNA
    elif valida(seq,RNA):
        return RNA
    elif valida(seq,PROT):
        return PROT
    else:
        raise TypeError("Sequencia invalida")
#%%

""" Funções que leem ficheiros e extraem a sequencia incluidas. A função ficheirotxt lé ficheiros
de texto e a função ficheirofasta lé ficheiros fasta
"""

def ficheirotxt(FileHandle):
    
    #verifica se a extensão do ficheiro é do tipo fasta
    x = re.search("txt$", FileHandle)
    
    if x:
        #caso for txt abre o ficheiro e lé o conteudo
        file = open(FileHandle, "r") 
        lines = file.readlines() 
        
        #string vazia a que vai ser adicionado as sequencias
        seq = ""
        for line in lines:
            #remove as quebras de linha no texto do ficheiro
            seq += line.replace("\n","") 
        file.close() 
        return seq
    
    #caso o input do ficheiro não for txt obtem-se um TypeError
    else:
        raise TypeError("Ficheiro invalido")
        

def ficheirofasta(FileHandle):
    
 
    #verifica se a extensão do ficheiro é do tipo fasta
    x = re.search("fasta$", FileHandle)
    
    
    if x:
        #caso for fasta abre o ficheiro e lé o conteudo
        ficheiro = open(FileHandle, "r") 
        ficheiro_fasta = ficheiro.read() 
        
        #Cria uma lista em que cada elemento é uma linha do ficheiro
        fasta_list = ficheiro_fasta.split("\n") 
        
        #Cria uma nova lista sem o primeiro elemnto da lista anterior (o cabeçalho do ficheiro fasta)
        seq_list = fasta_list[1:] 
        
        #Junta os elemntos da lista anterior para criar a sequencia
        seq = "".join(seq_list) 
        return seq
    
    #caso o input do ficheiro não for fasta obtem-se um TypeError
    else:
        raise TypeError("Ficheiro invalido")


#%%

"""
Funções que realizam diversas analises ás sequencias.

"""

def complemento_inverso(seq):
    
    #Verfica se é uma sequencia valida de DNA
    if valida(seq, DNA):
        
        #Inverte a sequencia e troca as bases 
        return (seq[::-1].upper().replace("A", "t").replace("T","a").replace("C","g").replace("G","c").upper())
    
    #Caso não for uma sequencia válida dá um TypeError
    else:
        raise TypeError("Sequencia invalida")

def transcricao(seq):
    
    #Verfica se é uma sequencia valida de DNA
    if valida(seq, DNA):
        
        #Troca as bases
        return seq.upper().replace("T", "U")
    
    #Caso não for uma sequencia válida dá um TypeError
    else:
        raise TypeError("Sequencia invalida")



def traducao(seq):
    
    #Verfica se é uma sequencia valida de DNA
    if valida(seq,DNA):
        
        #Dicionario com codões e respetivos aminoacidos
        gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        
        
        #Divide a sequencia em codoes
        seq_correct = re.findall('...',seq.upper())
        
        #string vazia que vai conter o resultado da tradução
        resultado = ""
        
        #Para cada codão na sequencia verifica se existe no dicionário e posteriormente adiciona o valor da chave na string vazia
        for codao in seq_correct:
            if codao in gencode.keys():
                resultado += gencode[codao]
        
        return resultado
    
    #Caso não for uma sequencia válida dá um TypeError
    else:
        raise TypeError("Sequencia invalida")
        
    
def contar_bases(seq):
    
    #dicionario vazio onde vai ser adicionado o valor da contagem das bases
    dic = {}
    
    #Invoca a função alfab para verificar que tipo de sequencia vai fazer a contagem
    if alfab(seq):
        
        #Faz a contagem de cada base existente na sequencia e adiciona ao dicionario
        for letter in seq.upper():
            if letter not in dic:
                dic[letter] = 1
            else:    
                dic[letter] +=1
                
        #Caso não haver certas bases é adicionado as bases que não foram encontradas ao dicionario e colocado o valor zero
        #Invoca a função alfab para verificar os caracteres permitidos no tipo de sequencia que foi verificada anteriomente
        for letter2 in alfab(seq):
            if letter2 not in dic:
                dic[letter2] = 0
                
        return dic
    
    #Caso não for uma sequencia válida dá um TypeError    
    else:
        raise TypeError("Sequencia invalida")
            



def reading_frames(seq):
    #lista vazia onde vai ser adicionado as reading frames
    lista = []
    
    #Verfica se é uma sequencia valida de DNA
    if valida(seq,DNA):
        
        #Adiciona 3 sequências traduzidas apartir do primeiro, segundo e terceiro nucleotido da sequência
        for i in range(3):
            lista.append(''.join(traducao(seq[i:])))
            
        #Adiciona a traducao do complemento inverso de 3 sequências apartir do primeiro, segundo e terceiro nucleotido da sequência original  
        for i in range (3):
            lista.append(traducao(complemento_inverso(seq)[i:]))
        
        return lista
    
    #Caso não for uma sequencia válida dá um TypeError 
    else:
        raise TypeError("Sequencia invalida")
            



def get_proteinas(seq):
    
    #Verfica se é uma sequencia valida de DNA
    if valida(seq,DNA):
        
        #Invoca a função anterior para ver as reading frames
        reading_frames(seq)
        
        #Extrai as proteinas nas reading frames procurando codões de inicio (o aminoacido M) até um codão stop ("_")
        prots = [re.findall('M.*?_', orf) for orf in reading_frames(seq)]
        
        #Cria uma lista que organiza as proteinas obtidas por tamanho e por ordem alfabética para as do mesmo tamanho
        lista = sorted(set(p for lista_p in prots for p in lista_p), key = lambda x: (-len(x), x))
        
        return lista
    #Caso não for uma sequencia válida dá um TypeError 
    else:
        raise TypeError("Sequencia invalida")








