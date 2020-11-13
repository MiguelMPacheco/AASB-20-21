# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:24:13 2020

@author: Dias
"""

from ficha4  import *
import unittest
import math


class fichaTest(unittest.TestCase):
    
    
    def teste_validar(self):
        
        
        #Teste de uma sequencia de DNA
        self.assertTrue(valida("ATGCTctcaTCGT", DNA), True) 
        
        #Teste de uma sequencia de RNA
        self.assertTrue(valida("AAAUUGccgcgc", RNA), True)
        
        #Teste de uma sequencia de proteina
        self.assertTrue(valida("LNRGGA_C", PROT), True)
        
        #Teste de uma sequencia invalida
        self.assertFalse(valida("ATATAxxxf", DNA), False)
        
    def teste_alfab(self):
        
        #Teste de uma sequencia de DNA
        self.assertEqual(alfab("ATATCGCGCG"), DNA)
        
        #Teste de uma sequencia de RNA
        self.assertEqual(alfab("AAAAUUUGCGCGC"), RNA)
        
        #Teste de uma sequencia de proteina
        self.assertEqual(alfab("LNRGGA_C"), PROT)
        
        #Teste de uma sequencia invalida
        with self.assertRaises(TypeError):
            resultado = alfab("ATACGUUULL")


    
    
    def teste_ficheiros(self):
        
        #Teste de leitura de ficheiros txt e fasta para uma sequencia de DNA
        self.assertEqual(ficheirotxt("sequencia_dna.txt"), "AAAATTACTCGATGAGATCGCTAGCGTGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGTTTAACCCTGTGTCAAGATCGCTAGCGTTGGTGTGTGGTGTGTGCTGCTGCAGATCGCTAGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGTTTAACCCTGTGTCAAGATCGGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGTTTAACCCTGTGTCAAGATCGCTAGCGTTGGTGTGCTAGCGTTGGTGTGGAGATCGCTAGCGTTGGTGTGAGATCGCTAGCGTTGGTGTGTTTAACCCTGTGTCAAGATCGCTAGCGTTGGTGTG")
        self.assertEqual(ficheirofasta("sequencia_dna.fasta"), "CTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTA")
        
        #Teste de leitura de ficheiros txt e fasta para uma sequencia de RNA
        self.assertEqual(ficheirotxt("sequencia_rna.txt"), "ACGCGCGACGUUCGACGCUGACGUUUCUACGUUCGUCUGACGUUGACGUUCGUCUGACGUUCUGGACGUUCGUCACGUUCGUCUUCGUCUGACCUGACGUUUCUACGCUGACGUGACGUUUUCUAGACGUUCGUUCGUCUGACGUUGACGUUCGUCUGACGUUCUGGACGUUCGUCACGUUCGUCUUCGUACGCGACGUUUGACGUUUCUACGUUCGUGACGUUCUGACGUUGACGUUCGUCUGACGUUCUGGACGUUCGUCACGUUCGUCGACGUUUUCGUCUGACGUUGACGUUCUGACGUUGACGUUCGUCUGACGUUCUGGACGUUCGUCACGUUCGUCGUUUCUACGUUCGUCUGACGUUGACGUUCGUCUGACGUUCUGGACGUUCGUCACGUUCGUC")
        self.assertEqual(ficheirofasta("sequencia_rna.fasta"), "AGACGCUUAUAUCUCGCUUAUAUCUCUCUGAUACGCUUAUAUUCUGAUACGCUACGGCUUAUAUCUCUCUGAUACGCUUAUAUACGCUUAUAUGUCUGAUACGCUCGCGAUUAUAUCUCUCUCUUAUAUGCGCGAUUAUAUCUCUCUUCUGAUACGCUUAUAUACGCUUAUAUGCGCGAUUAUAUGCUUAUAUCUCUCUGAUACGCUUAUAUACGCUUAUAUGCGCGAUUAUAUCUCUCUCUCUCUGAU")
        
        #Teste de leitura de ficheiros txt e fasta para uma sequencia de aminoacidos
        self.assertEqual(ficheirotxt("sequencia_prota.txt"), "YPWTQRFYPWRFFRFFDGGETYPWTTLGRLLVVGHFTEEDKATITSLWGKYPWTGETLGRLLVVYPWTETLGRLGETLGRLLVVYPWTRLLVVYPWTGETLGRLLVVYPWTGETLGVYPWTGETLGRLLVVGHFTEEDKATITSLWGKYPWTETLGRLGETLGRLLVVYPWTRLLVVYPWTGETLGRLLVVYPWTGETLGRYPWTTLLVVYPWTTQRFFDGGETYPWTTLGRLLVVYPWTGETLGRLLVVYPWTETLGRLGETLGRLLVVYPWTRLLVVYPWTGETLGRLLVVYPWTGETRLLDGGETLGRLLVVGHFTEEDKATITSLWGKGHFTEEDKATITSLWGKYPWTGETLGRLLVVYPWTETLGRLGETLGRLLVVYPWTRLLVVYPWTGETLGRLLVVYPWTGETLGRYPWTTLLVVYPWTTQRFFDGGETYPWTTLGRLLVVYPWTGETLGRLLVVYPWTETLGRLGETLGRLLVVYPWTRLLVVYPWTGETLGRLLVVYPWTGETLGRLLVVYPWT")
        self.assertEqual(ficheirofasta("sequencia_prota.fasta"), "MDSSEVVKVKQASIPAPGSILSQPNTEQSPAIVLPFQFEATTFGTAETAAQVSLQTADPITKLTAPYRHAQIVECKAILTPTDLAVSNPLTVYLAWVPANSPATPTQILRVYGGQSFVLGGAISAAKTIEVPLNLDSVNRMLKDSVTYTDTPKLLAYSRAPTNPSKIPTASIQISGRIRLSKPMLIAN")
        
        
        #Testes para verificar se a função identifica corretamente o tipo de ficheiro
        #Caso o ficheiro incorreto for colocado acontece um TypeError
        
        with self.assertRaises(TypeError):
            resultado = ficheirotxt("sequencia_prota.fasta")
            
        with self.assertRaises(TypeError):
            resultado = ficheirofasta("sequencia_prota.txt")
   
    
    def teste_complemento(self):
       
        #Teste de uma sequencia de DNA 
        self.assertEqual(complemento_inverso("ATGGTTTCA"), "TGAAACCAT") 
        
        #Testa se a sequência apresentar nucleotidos em minuscula
        self.assertEqual(complemento_inverso("ATGGtttCA"), "TGAAACCAT") 
        
        #Não faz complemento inverso caso não for uma sequencia valida de DNA
        with self.assertRaises(TypeError):
            resultado = complemento_inverso("uuugta")

    
    def teste_transcricao(self):
        
        
        #Teste de uma sequencia de DNA 
        self.assertEqual(transcricao("ATGCTCTCATCGT"), "AUGCUCUCAUCGU")
        
        # testa se a sequência apresentar nucleotidos em minuscula
        self.assertEqual(transcricao("ATGCTctcaTCGT"), "AUGCUCUCAUCGU")
        self.assertEqual(transcricao('atgc'), 'AUGC')
        
        #Se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
            resultado1 =  transcricao("ATttcgxxx")
        
        #Se a sequência já for de rna
        with self.assertRaises(TypeError):    
            resultado2 = transcricao('UAGUCUGUACuG') 

   
            
    
    def teste_traducao(self):
        
        #Teste de uma sequencia de DNA 
        self.assertEqual(traducao("TTAAATCGCGGTGGAGCATGTG"), "LNRGGAC") 
        
        #Testa se a sequência apresentar nucleotidos em minuscula
        self.assertEqual(traducao("TTAAATCGCGGTggaGCATGTG"), "LNRGGAC") 
        
        #Teste de uma sequencia de DNA caso a sequencia não for multipla de 3
        self.assertEqual(traducao("ATATAGC"), "I_")
        self.assertEqual(traducao("ATATAGCC"), "I_")
        self.assertEqual(traducao("AT"), "")
        
        #Se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
                resultado1 =  traducao("tatatafxxx")
                
    
    
    def teste_contar_bases(self):
        
        #Teste de uma sequencia de DNA
        self.assertEqual(contar_bases("ATATACGCGCG"), {'A': 3, 'T': 2, 'C': 3, 'G': 3})
        
        #Teste de uma sequencia de RNA
        self.assertEqual(contar_bases("AuAuAuCGCGCG"), {'A': 3, 'U': 3, 'C': 3, 'G': 3})
        
        #Teste de uma sequencia de aminoacidos
        self.assertEqual(contar_bases("LNRGGA_C"), {'L': 1, 'N': 1, 'R': 1,
                                                    'G': 2,'A': 1,'_': 1,'C': 1,'D': 0,'E': 0,
                                                    'F': 0,'H': 0,'I': 0,'K': 0,'M': 0,
                                                    'P': 0,'Q': 0,'S': 0,'T': 0,'V': 0,'W': 0,'Y': 0})
        
        
        #Teste de uma sequencia de DNA que não apresenta um ou mais tipos de bases(Tambem testa se a sequencia for em minuscula)
        self.assertEqual(contar_bases("atatat"), {'A': 3, 'T': 3, 'C': 0, 'G': 0} )

        #Se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
                resultado =  contar_bases("tatatafxxx")
    
    def teste_reading_frames(self):
        #Se apresentar sequência que não seja de DNA
        self.assertRaises(TypeError, reading_frames('XXXWWWJJKLMANTA'))
        self.assertRaises(TypeError, reading_frames('AUUUUTCCGGCGUAU'))
        self.assertRaises(TypeError, reading_frames('MGACDKWWTSY_'))
        #Se a sequência apresentar letras minúsculas
        self.assertEqual(reading_frames('atgCTGCATAtCTTTTAGCAagtGTCAGTAATAG'), ['MLHIF_QVSVI', 'CCISFSKCQ__', 'AAYLLASVSN', 'LLLTLAKRYAA', 'YY_HLLKDMQH', 'ITDTC_KICS']
                         
        self.assertEqual(reading_frames('ATGCTGCATATTCTTTTAGCAATGTCAGTCAACTAG'), ['MLHILLAMSVN_', 'CCIFF_QCQST', 'AAYSFSNVSQL', 'LVD_HC_KNMQH', '_LTDIAKRICS', 'S_LTLLKEYAA']
        pass
    
    def get_proteinas(self):
        #Caso a sequência seja inválida/ RNA/ Aminoácidos
        self.assertRaises(TypeError, get_proteins('AXKWYJZ'))
        self.assertRaises(TypeError, get_proteins('AUCGUUUA'))
        self.assertRaises(TypeError, get_proteins('KLMQP_'))
        #Testa letras pequenas e grandes
        self.assertEqual(get_proteins("atgCTGCATAtCTTTTAGCAagtGTCAGTAATAG"), ['MLHIF_'])
        #Caso a sequência tenha múltiplas proteínas
        self.assertEqual(get_proteins("ATGTCTGCGAATTAAGGGATGGTTTGTGAACGGTATAAATGAGCGCGCATGTTTAAGATACACGGGAACTAGGCGCGCTATAAATAATTTATGTAA"), ['MRAHLYRSQTIP_', 'MFKIHGN_','MVCERYK_','MSAHV_','MSAN_','M_'] )
        #Caso não tenha codões START/STOP
        self.assertRaises(TypeError, get_proteins('AGTCTCGCGCGTTTTTTAAATCGCGTTTTTTAG'))
        self.assertRaises(TypeError, get_proteins('ATGCTCGCGCGTTTTTTAAATCGCGTTTTTTAT'))
        pass
        
        
        




if __name__ == '__main__':
    unittest.main()
    
