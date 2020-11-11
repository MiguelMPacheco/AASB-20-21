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
        self.assertTrue(valida("ATGCTctcaTCGT"), True) 
        
        #Teste de uma sequencia de RNA
        self.assertFalse(valida("ATATAUUG"), False)
        
        #Teste de uma sequencia de proteina
        self.assertFalse(valida("LNRGGA_C"), False)
        
        #Teste de uma sequencia invalida
        self.assertFalse(valida("ATATAxxxf"), False)
        

    
    def teste_ficheirotxt(self):
        pass
    
    def teste_ficheirofasta(self):
        pass
    
    
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
        
        # se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
            resultado1 =  transcricao("ATttcgxxx")
        
        # se a sequência já for de rna
        with self.assertRaises(TypeError):    
            resultado2 = transcricao('UAGUCUGUACuG') 

   
            
    def teste_codoes(self):
        
        # testa uma sequencia com um numero de aminoacidos multiplo de 3
        self.assertEqual(codoes("TTAAATCGCGGTGGAGCATGTG"), "TTA AAT CGC GGT GGA GCA TGT") 
        #testa se o input tem letras minusculas e as converte em maiuscula
        self.assertEqual(codoes("TTAAATCGCGGTGGAGCATgtg"), "TTA AAT CGC GGT GGA GCA TGT") 
        #testa os casos em que o numero de codoes não é multiplo de 3
        self.assertEqual(codoes("ATATAGC"), "ATA TAG")
        self.assertEqual(codoes("ATATAGCC"), "ATA TAG")
        self.assertEqual(codoes("AT"), "")
        
        # se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
                resultado1 =  codoes("tatatafxxx")
    
    def teste_traducao(self):
        
        #Teste de uma sequencia de DNA 
        self.assertEqual(traducao("TTAAATCGCGGTGGAGCATGTG"), "LNRGGAC") 
        #Testa se a sequência apresentar nucleotidos em minuscula
        self.assertEqual(traducao("TTAAATCGCGGTggaGCATGTG"), "LNRGGAC") 
        
        #Teste de uma sequencia de DNA caso a sequencia não for multipla de 3
        self.assertEqual(traducao("ATATAGC"), "I_")
        self.assertEqual(traducao("ATATAGCC"), "I_")
        self.assertEqual(traducao("AT"), "")
        
        # se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
                resultado1 =  traducao("tatatafxxx")
                
    
    
    def teste_contar_bases(self):
        
        #Teste de uma sequencia de DNA
        self.assertEqual(contar_bases("ATATACGCGCG"), {'A': 3, 'T': 2, 'C': 3, 'G': 3})
        
        #Teste de uma sequencia de DNA que não apresenta um ou mais tipos da bases(Tambem testa se a sequencia for em minuscula)
        self.assertEqual(contar_bases("atatat"), {'A': 3, 'T': 3, 'C': 0, 'G': 0} )
        
        #Se a sequência apresentar elementos que não sejam nucleotidos
        with self.assertRaises(TypeError):
                resultado =  contar_bases("tatatafxxx")
    
    def teste_reading_frames(self):
        pass
    
    def get_proteinas(self):
        pass
        
        
        




if __name__ == '__main__':
    unittest.main()
    