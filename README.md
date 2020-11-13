# Resolução da ficha 4 de Algoritmos para Análise de Sequências Biológicas 20/21.

Esta resolução consiste em 2 ficheiros python assim como de 6 ficheiros de texto.
O ficheiro ficha4.py apresenta as seguintes funções para analisar sequências. Além das funções pedidas foi adicionado a função def alfab e a função valida foi alterada para validar qualquer tipo de sequência. 

### valida(seq, alf) 
Função que verifica se uma sequência é válida. Tens um argumento adicional alf em que se coloca o tipo de sequência que se quer verificar (DNA, RNA ou PROT). Devolve True ou False dependendo da sequencia

### def alfab(seq)
Função que devolve o tipo de sequência que se está a analisar (DNA, RNA, PROT)

### ficheirotxt(FileHandle) 
Função que recebe um ficheiro de texto aberto e que devolve uma sequência

### ficheirofasta(FileHandle) 
Função que recebe um ficheiro fasta aberto e que devolve uma sequência

### complemento_inverso(seq) 
Função que devolve o complemento inverso de uma sequência de DNA

### transcricao(seq) 
Função que devolve a transcrição de uma sequência de DNA

### traducao(seq) 
Função que devolve a tradução de uma sequência de DNA

### contar_bases(seq)
Função que conta as bases de uma sequência e devolve um dicionário com a contagem

### reading_frames(seq)
Função que devolve uma lista com as reading frames

### get_proteins(seq) 
Função que devolve a lista de todas as proteínas ordenadas por tamanho e por ordem alfabética para as do mesmo tamanho

O ficheiro testes_ficha4 realiza exaustivamente todo o tipo de testes nas funções apresentadas anteriormente para comprovar a veracidade destas. Foram anexadas neste trabalho 3 ficheiros txt e 3 ficheiros fasta retirados do NCBI para comprovar as funções de leitura de ficheiros.
