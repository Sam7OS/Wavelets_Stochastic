import numpy as np
import requests
import matplotlib.pyplot as plt

#Cette fonction permet de mettre un nombre renseigné de minimums à 0
def min_recursive_signal(signal,acc,nb_min):
    if(acc==nb_min):
        signal[signal==100000000]=0
        return signal
    else:
        signal[(np.abs(signal)==np.min(np.abs(signal)))]=100000000
        return min_recursive_signal(signal,(acc+1),nb_min)

def division_tab(signal,nb_periods,acc1,acc2):
    acc2=np.append(acc2,signal[acc1:(acc1+nb_periods)])
    if(acc1==len(signal)-nb_periods):
        return acc2
    else:
        return division_tab(signal,nb_periods,(acc1+1),acc2)

######################################################################################################################
#Transformée en ondelette de daubechie
######################################################################################################################

#La fonction daubechie_transform_recursive permet de calculer la transformée en ondelette de daubechie d'un signal
#Sa particularité est sa récursivité
def daubechie_transform_recursive(signal,acc1,acc2,acc3):
    #Il s'agit des coefficients de daubechie à 4 points
    daub_coef_somme=np.array([0.4829629,0.8365163,0.2241439,-0.1294095])
    daub_coef_soustraction=np.array([-0.1294095,-0.2241439,0.8365163,-0.4829629])
    #Nous séparons les cas où la longueur du signal est supérieur à 2 et le cas où elle vaut 2 
    if(len(signal)>2):
        #Nous rajoutons à acc2 le résultat du produit entre 4 points du signal et les coefficients de somme
        acc2=np.append(acc2,np.dot(signal[acc1:(acc1+len(daub_coef_somme))],daub_coef_somme))
        #Nous rajoutons à acc3 le résultat du produit entre 4 points du signal et les coefficients de soustraction
        acc3=np.append(acc3,np.dot(signal[acc1:(acc1+len(daub_coef_somme))],daub_coef_soustraction))
        """
        Nous incrémentons acc1 jusqu'à ce qu'il soit égal à la différence 
        entre la longueur du signal et le nombre de coefficients de daubechie
        """
        if(acc1==(len(signal)-len(daub_coef_somme))):
            """
            Si acc1 est égal à la longueur du signal moins le nombre de coefficients
            alors nous prenons les deux derniers et premiers éléments du signal afin
            de les concaténer
            Puis nous rajoutons à acc2 et acc3 le résultat du produit entre ce tableau concaténé
            et les coefficients de daubechie pour la somme et la différence
            Enfin nous retournons la concaténation de acc2 et acc3 car si la condition de cette boucle est 
            vérifiée, cela veut dire que nous avons parcouru tout le signal
            """
            z=np.concatenate((signal[-(len(daub_coef_somme)//2):len(signal)],signal[0:(len(daub_coef_somme)//2)]))
            acc2=np.append(acc2,np.dot(z,daub_coef_somme))
            acc3=np.append(acc3,np.dot(z,daub_coef_soustraction))
            return np.concatenate((acc2,acc3))
        else:
            """
            Si la condition n'est pas vérifiée alors cela veut dire que nous ne sommes pas arrivés
            au bout de la liste et donc il faut continuer, c'est pourquoi nous rappelons de manière récursive 
            la fonction afin de répéter les mêmes opérations sur le segment suivant du signal
            """
            return daubechie_transform_recursive(signal,(acc1+2),acc2,acc3)
    else:
        """
        Si la longueur du signal vaut 2 alors nous sommes dans un cas particulier:
        Nous devons dans ce cas, concaténer le signal avec lui-même afin qu'il ait une longueur égale au 
        nombre de coefficients de daubechie puis nous rajoutons à nouveau le résultat du produit entre
        ce tableau concaténé et les coefficients de daubechie pour la somme et la différence
        Enfin nous retournons la concaténation de acc2 et acc3 car nous avons déjà parcouru tout le signal 
        """
        z2=np.concatenate((signal,signal))
        acc2=np.append(acc2,np.dot(z2,daub_coef_somme))
        acc3=np.append(acc3,np.dot(z2,daub_coef_soustraction))
        return np.concatenate((acc2,acc3))

#La fonction daubechie_transform permet de calculer la transformée en ondelette de daubechie d'un signal
def daubechie_transform(signal):
    """
    Cette fonction sert juste à fixer les paramètres de la fonction daubechie_transform_recursive
    afin que son utilisation par la suite soit plus aisé. De ce fait, pour faire une transformée en 
    ondelette, nous n'aurons qu'à appeler cette fonction avec le signal renseigné en paramètre
    """
    indice=0
    tab_somme=np.array([])
    tab_soustraction=np.array([])
    return daubechie_transform_recursive(signal,indice,tab_somme,tab_soustraction)

"""
La fonction daubechie_transform_inverse_recursive permet de calculer la matrice de daubechie 
nécessaire à la transformée inverse d'un signal.Tout comme la fonction de transformée, elle est récursive.
"""
def daubechie_transform_inverse_recursive(z,acc1,acc1b,acc2,acc3):
    daub_coef_somme=np.array([0.4829629,0.8365163,0.2241439,-0.1294095])
    daub_coef_soustraction=np.array([-0.1294095,-0.2241439,0.8365163,-0.4829629])
    #Nous séparons les cas où acc1b est inférieur à la longueur du signal divisé par 2 et le cas où il y est égal
    if(acc1b<(len(z)//2)):
        """
        Nous vérifions d'abord que acc1 est inférieur à la longueur du signal 
        moins le nombre de coefficients divisé par 2.
        """
        if(acc1==(len(z)-(len(daub_coef_somme)//2))):
            """
            Si la condition est vérifié alors nous sommes arrivés à l'avant-dernière ligne de la matrice
            et donc il nous faut remplir la dernière ligne pour laquelle les 2 dernières colonnes correspondent
            aux deux premiers coefficients de daubechie et les deux premières colonnes correspondent aux deux
            derniers coefficients de daubechie. Nous rappelons par la suite la fonction afin que puisse être
            vérifiée la première condition et donc retourner la concaténation de acc2 et acc3
            """
            acc2[acc1b][(len(acc2[acc1b])-(len(daub_coef_somme)//2)):(len(acc2[acc1b]))]=daub_coef_somme[0:(len(daub_coef_somme)//2)]
            acc2[acc1b][0:(len(daub_coef_somme)//2)]=daub_coef_somme[(len(daub_coef_somme)//2):(len(daub_coef_somme))]            
            acc3[acc1b][(len(acc3[acc1b])-(len(daub_coef_somme)//2)):(len(acc3[acc1b]))]=daub_coef_soustraction[0:(len(daub_coef_soustraction)//2)]
            acc3[acc1b][0:(len(daub_coef_somme)//2)]=daub_coef_soustraction[(len(daub_coef_somme)//2):(len(daub_coef_somme))]
            return daubechie_transform_inverse_recursive(z,(acc1+2),(acc1b+1),acc2,acc3)
        else:
            """
            Si nous ne sommes pas arrivés à l'avant-dernière ligne de la matrice alors nous continuons
            de remplir la matrice avec les coefficients de daubechie. La matrice étant à la base remplie 
            de zéro, nous remplissons chaque ligne avec les 4 coefficients décalés de 2 colonnes à chaque fois.
            Nous rappelons ensuite la fonction afin d'appliquer les mêmes étapes à la ligne suivante de la
            matrice.
            """
            acc2[acc1b][acc1:(acc1+len(daub_coef_somme))]=daub_coef_somme
            acc3[acc1b][acc1:(acc1+len(daub_coef_somme))]=daub_coef_soustraction
            return daubechie_transform_inverse_recursive(z,(acc1+2),(acc1b+1),acc2,acc3)
    else:
        #Lorsque acc1b vaut enfin la longueur du signal divisé par 2 alors nous pouvons retourner la matrice de daubechie
        return np.concatenate((acc2,acc3),axis=0)

#Cette fonction permet de calculer la transformée en ondelette inverse
def daubechie_transform_inverse(z):
    result=np.array([])
    if(len(z)>2):
        """
        Si nous sommes dans un cas où la longueur du signal est supérieur à 2 alors nous 
        fixons les paramètres de la fonction daubechie_transform_inverse_recursive et
        nous calculons la matrice de daubechie adéquate au signal. Enfin, nous retournons le
        produit matriciel entre le signal et la matrice de daubechie.
        """
        A=np.zeros(((len(z)//2),len(z)))
        B=np.zeros(((len(z)//2),len(z)))
        indice_ligne=0
        indice_colonne=0
        daubechie_matrix=daubechie_transform_inverse_recursive(z,indice_ligne,indice_colonne,A,B)
        result=np.dot(z,daubechie_matrix)
    else:
        """
        Si la longueur du signal est égale à 2 nous concaténons le signal avec lui-même puis
        nous répétons les mêmes opérations que dans la boucle précédente sur ce signal concaténé
        """
        z2=np.concatenate((z,z,z))
        A=np.zeros(((len(z2)//2),len(z2)))
        B=np.zeros(((len(z2)//2),len(z2)))
        indice_ligne=0
        indice_colonne=0
        daubechie_matrix=daubechie_transform_inverse_recursive(z2,indice_ligne,indice_colonne,A,B)
        result=np.dot(z2,daubechie_matrix)[0:2]
    return result

#Cette fonction permet de débruiter un signal
def smooth_daubechie(signal,pourcentage_nb):
    #Tout d'abord, nous calculons la transformée en ondelette du signal
    z=daubechie_transform(signal)
    """
    Puis nous mettons à 0 tous les éléments du signal transformée dont la valeur absolue 
    est inférieur à un seuil arbitraire
    """
    z_deb=min_recursive_signal(z,0,np.ceil(pourcentage_nb*len(signal)/100))
    #Enfin nous retournons le signal débruité après avoir fait une transformée inverse
    return daubechie_transform_inverse(z_deb)

#Cette fonction permet de faire une transformée en ondelette multiple. Elle est aussi récursive.
def two_level_daubechie_transform(signal,level,acc):
    #Nous commençons d'abord par faire une première transformée du signal
    z=daubechie_transform(signal)
    #Puis nous rajoutons le résultat de cette transformée à l'accumulateur
    acc=np.append(acc,z)
    #Nous regardons si nous avons atteint le niveau maximal
    if(level!=0):
        """
        Si nous n'avons pas atteint le niveau maximal alors nous rappelons la fonction sur la première moitié
        du signal et ainsi de suite. Cela permet de rajouter à l'accumulateur tous les résultats de ces transformées
        """
        return two_level_daubechie_transform(np.split(z,2)[0],(level-1),acc)
    else:
        """
        Lorsque nous avons atteint le niveau maximal alors nous renvoyons acc qui contient les résultats 
        de toutes les transformées successives
        """
        return acc

#Cette fonction permet de trier les résultats obtenus avec la fonction two_level_daubechie_transform
def multi_daubechie(signal,level):
    #Tout d'abord, nous appelons la fonction two_level_daubechie_transform
    z=two_level_daubechie_transform(signal,level,np.array([]))
    #Puis nous calculons des puissances de 2 
    powers=np.arange((level+1),0,-1)
    puiss_2=np.power(2,powers)
    #Nous faisons une somme cumulative de ces puissances de 2
    indexs=np.cumsum(puiss_2)
    #Enfin nous divisons le tableau z aux indices contenus dans le tableau indexs
    z_s=np.split(z,indexs)
    #Nous retournons le tableau splitté en plusieurs sous-tableaux
    return z_s

"""
Cette fonction permet de récupérer les éléments du signal transformé plusieurs fois 
grâce au tableau renvoyé par la fonction multi_daubechie
"""
def multi_trans_daubechie(z,level):
    """
    Cette fonction permet de sélectionner les éléments
    qui correspondent à la transformée multiple du signal
    """
    z[(level-1)][0:(len(z[level]))]=z[level]
    if(level!=1):
        return multi_trans_daubechie(z,(level-1))
    else:
        return z[0]

#Cette fonction permet de faire une transformée multiple à l'aide des fonctions écrites au-dessus
def transformee_daubechie_multi(signal,level):
    #Dans un premier temps, nous calculons toutes les transformées possibles
    z=multi_daubechie(signal,level)
    #Puis nous retournons uniquement la transformée multiple au niveau le plus élevé
    return multi_trans_daubechie(z,level)

#Cette fonction permet de splitter le signal transformé plusieurs fois
def trans_daubechie_multi_inverse(signal,level):
    #Nous calculons les puissances de 2
    powers=np.arange(1,(level+2),1)
    puiss_2=np.power(2,powers)
    #Puis nous splittons le signal transformé avec comme indices les puissances de 2
    z_s=np.split(signal,puiss_2)
    #Nous retournons enfin le tableau du signal splitté en sous-tableau
    return z_s

#Cette fonction permet de faire une transformée multiple inverse et donc de revenir au signal d'origine
def trans_daubechie_inverse(signal,level,acc,acc2):
    #Tout d'abord, nous splittons le signal transformé
    z=trans_daubechie_multi_inverse(signal,level)
    #Ensuite, nous rajoutons à acc le bon sous-tableau du signal splitté
    acc=np.append(acc,z[acc2])
    #Enfin, nous appliquons une transformée inverse à acc
    acc=daubechie_transform_inverse(acc)
    if(level!=acc2):
        #Si nous n'avons pas atteint le niveau maximal alors nous devons continuer et nous rappelons la fonction
        return trans_daubechie_inverse(signal,level,acc,(acc2+1))
    else:
        #Si nous avons atteint le niveau maximal, alors nous renvoyons le signal d'origine
        return acc

#Cette fonction permet de débruiter un signal à l'aide d'une transformée en ondelette multiple
def smooth_daubechie_multi(signal,level,pourcentage_nb):
    #Tout d'abord, nous calculons la transformée en ondelette multiple du signal
    z=transformee_daubechie_multi(signal,level)
    """
    Puis nous mettons à 0 tous les éléments du signal transformée dont la valeur absolue 
    est inférieur à un seuil arbitraire
    """
    z_deb=min_recursive_signal(z,0,np.ceil(pourcentage_nb*len(signal)/100))
    #Enfin nous retournons le signal débruité après avoir fait une transformée inverse
    return trans_daubechie_inverse(z_deb,level,np.array([]),0)

######################################################################################################################
#Transformée en ondelette de Haar
######################################################################################################################

#Cette fonction permet de faire une transformée avec une ondelette de Haar
def w_haar(signal):
    #Voici les coefficients de Haar
    haar_coef_somme=np.array([(1/np.sqrt(2)),(1/np.sqrt(2))])
    haar_coef_soustraction=np.array([(1/np.sqrt(2)),(-1/np.sqrt(2))])
    #Nous divisions le signal en plusieurs pairs de nombres
    pairs=np.split(signal,(len(signal)//(len(haar_coef_somme))))
    #Nous faisons un produit matriciel de ces pairs avec les coefficient de Haar
    tab_somme=np.dot(haar_coef_somme,np.transpose(pairs))
    tab_soustraction=np.dot(haar_coef_soustraction,np.transpose(pairs))
    #Nous renvoyons enfin le résultat de ces produits concaténés
    return np.concatenate((tab_somme,tab_soustraction))

#Cette fonction permet de faire une transformée inverse avec une ondelette de Haar  
def w_haar_inverse(signal):
    #Voici les coefficients de Haar
    haar_coef_somme=np.array([(1/np.sqrt(2)),(1/np.sqrt(2))])
    haar_coef_soustraction=np.array([(1/np.sqrt(2)),(-1/np.sqrt(2))])
    #Nous divisons le signal transformée en plusieurs singletons
    coefs=np.split(signal,(len(signal)))
    #Nous multiplions chaque singleton de la première moitié du tableau splitté par les coefficients de somme de Haar 
    tab_somme_inverse=coefs[0:(len(signal)//2)]*haar_coef_somme
    #Nous multiplions chaque singleton de la seconde moitié du tableau splitté par les coefficients de différence de Haar     
    tab_soustraction_inverse=coefs[(len(signal)//2):len(signal)]*haar_coef_soustraction
    #Nous additionnons les deux matrices obtenues
    result=(tab_soustraction_inverse+tab_somme_inverse)
    #Nous pouvons ensuite retourner le signal d'origine
    result=np.split(result,(len(signal)//2))
    return np.concatenate((result),axis=1)[0]

#Cette fonction permet de débruiter un signal
def smooth_haar(signal,pourcentage_nb):
    #Nous faisons d'abord une transformée du signal
    z=w_haar(signal)
    #Puis nous mettons à 0 les valeurs qui sont inférieurs à un certain seuil
    z_deb=min_recursive_signal(z,0,np.ceil(pourcentage_nb*len(signal)/100))
    #Enfin nous retournons la transformée inverse du signal débruité
    return w_haar_inverse(z_deb)

#Cette fonction permet de faire une transformée en ondelette multiple. Elle est récursive.
def two_level_haar_transform(signal,level,acc):
    #Nous commençons d'abord par faire une première transformée du signal
    z=w_haar(signal)
    #Puis nous rajoutons le résultat de cette transformée à l'accumulateur
    acc=np.append(acc,z)
    #Nous vérifions si nous avons atteint le niveau maximal
    if(level!=0):
        """
        Si nous n'avons pas atteint le niveau maximal alors nous rappelons la fonction sur la première moitié
        du signal et ainsi de suite. Cela permet de rajouter à l'accumulateur tous les résultats de ces transformées
        """
        return two_level_haar_transform(np.split(z,2)[0],(level-1),acc)
    else:
        """
        Lorsque nous avons atteint le niveau maximal alors nous renvoyons acc qui contient les résultats 
        de toutes les transformées successives
        """
        return acc

#Cette fonction permet de trier les résultats obtenus avec la fonction two_level_haar_transform
def multi(signal,level):
    #Tout d'abord, nous appelons la fonction two_level_haar_transform
    z=two_level_haar_transform(signal,level,np.array([]))
    #Puis nous calculons des puissances de 2 
    powers=np.arange((level+1),0,-1)
    puiss_2=np.power(2,powers)
    #Nous faisons une somme cummulative de ces puissances de 2
    indexs=np.cumsum(puiss_2)
    #Enfin nous divisons le tableau z aux indices contenus dans le tableau indexs
    z_s=np.split(z,indexs)
    #Nous retournons le tableau splitté en plusieurs sous-tableaux
    return z_s

"""
Cette fonction permet de récupérer les éléments du signal transformé plusieurs fois 
grâce au tableau renvoyé par la fonction multi
"""
def multi_trans(z,level):
    """
    Cette fonction permet de sélectionner les éléments
    qui correspondent à la transformée multiple du signal
    """
    z[(level-1)][0:(len(z[level]))]=z[level]
    if(level!=1):
        return multi_trans(z,(level-1))
    else:
        return z[0]

#Cette fonction permet de faire une transformée multiple à l'aide des fonctions écrites au-dessus
def transformee_multi(signal,level):
    #Dans un premier temps, nous calculons toutes les transformées possibles
    z=multi(signal,level)
    #Puis nous retournons uniquement la transformée multiple au niveau le plus élevé
    return multi_trans(z,level)

#Cette fonction permet de splitter le signal transformé plusieurs fois
def trans_multi_inverse(signal,level):
    #Nous calculons les puissances de 2
    powers=np.arange(1,(level+2),1)
    puiss_2=np.power(2,powers)
    #Puis nous splittons le signal transformé avec comme indices les puissances de 2
    z_s=np.split(signal,puiss_2)
    #Nous retournons enfin le tableau du signal splitté en sous-tableau
    return z_s

#Cette fonction permet de faire une transformée multiple inverse et donc de revenir au signal d'origine
def trans_inverse(signal,level,acc,acc2):
    #Tout d'abord, nous splittons le signal transformé
    z=trans_multi_inverse(signal,level)
    #Ensuite, nous rajoutons à acc le bon sous-tableau du signal splitté
    acc=np.append(acc,z[acc2])
    #Enfin, nous appliquons une transformée inverse à acc
    acc=w_haar_inverse(acc)
    if(level!=acc2):
        #Si nous n'avons pas atteint le niveau maximal alors nous devons continuer et nous rappelons la fonction
        return trans_inverse(signal,level,acc,(acc2+1))
    else:
        #Si nous avons atteint le niveau maximal, alors nous renvoyons le signal d'origine
        return acc

#Cette fonction permet de débruiter un signal à l'aide d'une transformée en ondelette multiple
def smooth_haar_multi(signal,level,pourcentage_nb):
    #Tout d'abord, nous calculons la transformée en ondelette multiple du signal
    z=transformee_multi(signal,level)
    """
    Puis nous mettons à 0 tous les éléments du signal transformée dont la valeur absolue 
    est inférieur à un seuil arbitraire
    """
    z_deb=min_recursive_signal(z,0,np.ceil(pourcentage_nb*len(signal)/100))
    #Enfin nous retournons le signal débruité après avoir fait une transformée inverse
    return trans_inverse(z_deb,level,np.array([]),0)

#Cette fonction permet de faire des prédictions avec une régression polynomiale
def prediction(x_train,y_train,x_test,degree):
    coefs=np.polyfit(x_train,y_train,degree)
    poly=np.poly1d(coefs)
    pred=poly(x_test)
    return pred
  
#####################################################################################################################
#La classe BTC_Price nous permet de récupérer les données du bitcoin directement sur Internet
#avec l'aide de l'API de Kraken, un exchange de devise en ligne 
#La paire de devise utilisé est EURO/Bitcoin
#####################################################################################################################

class BTC_Price:
    def __init__(self,time_interval):
        response=requests.get("https://api.kraken.com/0/public/OHLC?pair=BTCEUR&interval="+str(time_interval))  
        data=response.json()
        listprices_high=np.zeros(len(data["result"]["XXBTZEUR"]))
        listprices_low=np.zeros(len(data["result"]["XXBTZEUR"]))
        list_time=np.zeros(len(data["result"]["XXBTZEUR"]))
        for k in range(0,len(data["result"]["XXBTZEUR"])):
            list_time[k]=((float(data["result"]["XXBTZEUR"][k][0])-float(data["result"]["XXBTZEUR"][0][0]))/60) 
            listprices_high[k]=float(data["result"]["XXBTZEUR"][k][2])
            listprices_low[k]=float(data["result"]["XXBTZEUR"][k][3])
        
        self.list_prices_max=listprices_high
        self.list_prices_min=listprices_low
        self.list_times=list_time
    
    #Cette fonction permet de calculer une moyenne mobile simple
    def SMA(self,nb_periods,Max):
        list_prices=self.list_prices_min
        if Max==True:
            list_prices=self.list_prices_max
        z=division_tab(list_prices,nb_periods,0,np.array([]))
        index=np.arange(0,nb_periods*len(list_prices),nb_periods)
        z_s=np.split(z,index)
        z_m=np.average(z_s[1:(len(z_s)-nb_periods+1)],axis=1)
        return z_m
    
    #Cette fonction permet de calculer les écarts-types mobiles
    def Std_M(self,nb_periods,Max):
        list_prices=self.list_prices_min
        if Max==True:
            list_prices=self.list_prices_max
        z=division_tab(list_prices,nb_periods,0,np.array([]))
        index=np.arange(0,nb_periods*len(list_prices),nb_periods)
        z_s=np.split(z,index)
        z_m=np.std(z_s[1:(len(z_s)-nb_periods+1)],axis=1)
        return z_m    
    
    #Cette fonction permet de calculer une moyenne mobile exponentielle
    def EMA(self,nb_periods,Max):
        list_prices=self.list_prices_min
        if Max==True:
            list_prices=self.list_prices_max
        z=division_tab(list_prices,nb_periods,0,np.array([]))
        index=np.arange(0,nb_periods*len(list_prices),nb_periods)
        z_s=np.split(z,index)
        poids=np.arange(0,nb_periods,1)
        alpha=np.power(2,poids)
        z_m=np.average(z_s[1:(len(z_s)-nb_periods+1)],axis=1,weights=alpha)
        return z_m
    
    #Cette fonction permet de calculer les bandes de Bollinger
    def Bollinger(self,nb_periods,Max):
        SMA=self.SMA(nb_periods,Max)
        Std=self.Std_M(nb_periods,Max)
        boll_1=SMA+2*Std
        boll_2=SMA-2*Std
        return boll_1,boll_2
            
    #Fonction utilisé pour afficher les courbes et les prédictions
    def Affiche(self,pct_debruit_daubechie,pct_debtruit_haar,nb_period_SMA,degree):
        a=511
        b=530

        interval_pred=self.list_times[a:b]
        courbe_reel=self.list_prices_max[a:b]

        x_train=self.list_times[0:512]
        y_train=self.list_prices_max[0:512]
        
        Boll=self.Bollinger(nb_period_SMA,True)
        
        predict_1=prediction(x_train,y_train,interval_pred,degree)
        
        y_train_2=smooth_daubechie_multi(self.list_prices_max[0:512],8,pct_debruit_daubechie)
        predict_2=prediction(x_train,y_train_2,interval_pred,degree)
        
        y_train_3=smooth_haar_multi(self.list_prices_max[0:512],8,pct_debtruit_haar)
        print(y_train_3)
        predict_3=prediction(x_train,y_train_3,interval_pred,degree)
                    
        y_train_4=self.SMA(nb_period_SMA,True)[0:512]
        predict_4=prediction(x_train,y_train_4,interval_pred,degree)
        
        plt.figure()
        plt.subplot(2, 2, 1)
        plt.title("Corrélation prédictions ="+str(np.corrcoef(predict_1,courbe_reel)[0][1]))
        plt.plot(x_train,y_train,label="courbe train")
        plt.plot(interval_pred,courbe_reel,label="courbe réelle")
        plt.plot(interval_pred,predict_1,label="prediction direct")
        plt.plot(x_train,Boll[0][0:512],label="Bollinger supérieur")
        plt.plot(x_train,Boll[1][0:512],label="Bollinger inférieur")
        plt.xlabel("temps en secondes")
        plt.ylabel("Courbes")
        plt.legend(loc="upper left")
        
        plt.subplot(2, 2, 2)
        plt.title("Corrélation prédictions ="+str(np.corrcoef(predict_2,courbe_reel)[0][1]))
        plt.plot(x_train,y_train,label="courbe train")
        plt.plot(interval_pred,courbe_reel,label="courbe réelle")
        plt.plot(x_train,y_train_2,label="transformée daubechie")
        plt.plot(interval_pred,predict_2,label="prediction daubechie")
        plt.xlabel("temps en secondes")
        plt.ylabel("Courbes")
        plt.legend(loc="upper left")
        
        plt.subplot(2, 2, 3)
        plt.title("Corrélation prédictions ="+str(np.corrcoef(predict_3,courbe_reel)[0][1]))
        plt.plot(x_train,y_train,label="courbe train")
        plt.plot(interval_pred,courbe_reel,label="courbe réelle")
        plt.plot(x_train,y_train_3,label="transformée haar")
        plt.plot(interval_pred,predict_3,label="prediction haar")
        plt.xlabel("temps en secondes")
        plt.ylabel("Courbes")
        plt.legend(loc="upper left")

        plt.subplot(2, 2, 4)
        plt.title("Corrélation prédictions ="+str(np.corrcoef(predict_4,courbe_reel)[0][1]))
        plt.plot(x_train,y_train,label="courbe train")
        plt.plot(interval_pred,courbe_reel,label="courbe réelle")
        plt.plot(x_train,y_train_4,label="SMA")
        plt.plot(interval_pred,predict_4,label="prediction SMA")
        plt.xlabel("temps en secondes")
        plt.ylabel("Courbes")
        plt.legend(loc="upper left")
        
        plt.show()

#Applications :
BTC=BTC_Price(60)
BTC.Affiche(85,85,26,15)