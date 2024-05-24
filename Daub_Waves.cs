namespace WaveLibrary;

using ArrayCalc;

public class Daub_Waves // Transformée en ondelette de daubechie
{
    // Il s'agit des coefficients de daubechie à 4 points
    private readonly Vector daub_coef_somme = new([0.4829629,0.8365163,0.2241439,-0.1294095]);
    private readonly Vector daub_coef_soustraction = new([-0.1294095,-0.2241439,0.8365163,-0.4829629]);

    private readonly int daub_len;

    public Vector trans_signal;

    public Vector invtrans_signal;

    public Vector smooth_daub;

    public Vector smooth_daub_multi;
    
    public Daub_Waves(Vector signal, int pct_smooth = 0){
        Console.WriteLine("/////////////////////////////// Daubechie ////////////////////////////////////////////");
        daub_len = daub_coef_somme.VecSize();
        Console.WriteLine("Début signal");
        signal.VecDisp();
        Console.WriteLine("Fin signal");
        trans_signal = Daubechie_transform(signal);
        Console.WriteLine("Début signal transformé");
        trans_signal.VecDisp();
        Console.WriteLine("Fin signal transformé");
        invtrans_signal = Daubechie_transform_inverse(trans_signal);
        Console.WriteLine("Début signal transformé inverse");
        invtrans_signal.VecDisp();
        Console.WriteLine("Fin signal transformé inverse");
        smooth_daub = Smooth_daubechie(signal,pct_smooth);
        Console.WriteLine("Début signal débruité");
        smooth_daub.VecDisp();
        Console.WriteLine("Fin signal débruité");
        smooth_daub_multi = Smooth_daubechie_multi(signal,(int)Math.Log2(signal.VecSize()-1),pct_smooth);
        Console.WriteLine("Début signal débruité multi");
        smooth_daub_multi.VecDisp();
        Console.WriteLine("Fin signal débruité multi");
        Console.WriteLine("//////////////////////////////////////////////////////////////////////////////////////");
    }

    // La fonction daubechie_transform_recursive permet de calculer la transformée en ondelette de daubechie d'un signal
    // Sa particularité est sa récursivité
    private Vector Daubechie_transform_recursive(Vector signal, int acc1, Vector acc2, Vector acc3){
        // Nous séparons les cas où la longueur du signal est supérieur à 2 et le cas où elle vaut 2
        if(signal.VecSize() > 2){
            // Nous rajoutons à tab_somme le résultat du produit entre 4 points du signal et les coefficients de somme
            acc2 = acc2.VecApp(signal.SubVec(acc1,acc1+daub_len).VecDot(daub_coef_somme));
            // Nous rajoutons à tab_soustraction le résultat du produit entre 4 points du signal et les coefficients de soustraction
            acc3 = acc3.VecApp(signal.SubVec(acc1,acc1+daub_coef_soustraction.VecSize()).VecDot(daub_coef_soustraction));
            /*
                Nous incrémentons acc1 jusqu'à ce qu'il soit égal à la différence 
                entre la longueur du signal et le nombre de coefficients de daubechie
            */
            if(acc1 == signal.VecSize()-daub_len){
                /*
                    Si acc1 est égal à la longueur du signal moins le nombre de coefficients
                    alors nous prenons les deux derniers et premiers éléments du signal afin
                    de les concaténer
                    Puis nous rajoutons à acc2 et acc3 le résultat du produit entre ce tableau concaténé
                    et les coefficients de daubechie pour la somme et la différence
                    Enfin nous retournons la concaténation de acc2 et acc3 car si la condition de cette boucle est 
                    vérifiée, cela veut dire que nous avons parcouru tout le signal
                */
                Vector z = signal.SubVec(signal.VecSize()-daub_len/2,signal.VecSize()).VecCon(signal.SubVec(0,daub_len/2));
                acc2=acc2.VecApp(z.VecDot(daub_coef_somme));
                acc3=acc3.VecApp(z.VecDot(daub_coef_soustraction));
                return acc2.VecCon(acc3);

            }
            else{
                /*
                    Si la condition n'est pas vérifiée alors cela veut dire que nous ne sommes pas arrivés
                    au bout de la liste et donc il faut continuer, c'est pourquoi nous rappelons de manière récursive 
                    la fonction afin de répéter les mêmes opérations sur le segment suivant du signal
                */
                return Daubechie_transform_recursive(signal,acc1+2,acc2,acc3);
            }
        }
        else{
            /*
                Si la longueur du signal vaut 2 alors nous sommes dans un cas particulier:
                Nous devons dans ce cas, concaténer le signal avec lui-même afin qu'il ait une longueur égale au 
                nombre de coefficients de daubechie puis nous rajoutons à nouveau le résultat du produit entre
                ce tableau concaténé et les coefficients de daubechie pour la somme et la différence
                Enfin nous retournons la concaténation de acc2 et acc3 car nous avons déjà parcouru tout le signal 
            */
            Vector z2 = signal.VecCon(signal);
            acc2 = acc2.VecApp(z2.VecDot(daub_coef_somme));
            acc3 = acc3.VecApp(z2.VecDot(daub_coef_soustraction));
            return acc2.VecCon(acc3);
        }
    }

    // La fonction daubechie_transform permet de calculer la transformée en ondelette de daubechie d'un signal
    private Vector Daubechie_transform(Vector signal){
        /*
            Cette fonction sert juste à fixer les paramètres de la fonction daubechie_transform_recursive
            afin que son utilisation par la suite soit plus aisé. De ce fait, pour faire une transformée en 
            ondelette, nous n'aurons qu'à appeler cette fonction avec le signal renseigné en paramètre
        */
        int indice=0;
        Vector tab_somme=new([]);
        Vector tab_soustraction=new([]);
        return Daubechie_transform_recursive(signal,indice,tab_somme,tab_soustraction);
    }

    /*
        La fonction daubechie_transform_inverse_recursive permet de calculer la matrice de daubechie 
        nécessaire à la transformée inverse d'un signal.Tout comme la fonction de transformée, elle est récursive.
    */
    private Matrix Daubechie_transform_inverse_recursive(Vector z, int acc1, int acc1b, Matrix acc2, Matrix acc3){
        //Nous séparons les cas où acc1b est inférieur à la longueur du signal divisé par 2 et le cas où il y est égal
        if(acc1b < z.VecSize()/2){
            /*
                Nous vérifions d'abord que acc1 est inférieur à la longueur du signal 
                moins le nombre de coefficients divisé par 2.
            */
            if(acc1 == z.VecSize()-daub_len/2){
                /*
                    Si la condition est vérifié alors nous sommes arrivés à l'avant-dernière ligne de la matrice
                    et donc il nous faut remplir la dernière ligne pour laquelle les 2 dernières colonnes correspondent
                    aux deux premiers coefficients de daubechie et les deux premières colonnes correspondent aux deux
                    derniers coefficients de daubechie. Nous rappelons par la suite la fonction afin que puisse être
                    vérifiée la première condition et donc retourner la concaténation de acc2 et acc3
                */
                acc2.Change_values(acc1b,acc2.Nb_cols()-daub_len/2,daub_coef_somme.SubVec(0,daub_len/2));
                acc2.Change_values(acc1b,0,daub_coef_somme.SubVec(daub_len/2,daub_len));   
                acc3.Change_values(acc1b,acc3.Nb_cols()-daub_len/2,daub_coef_soustraction.SubVec(0,daub_len/2));
                acc3.Change_values(acc1b,0,daub_coef_soustraction.SubVec(daub_len/2,daub_len));
                return Daubechie_transform_inverse_recursive(z,acc1+2,acc1b+1,acc2,acc3);
            }
            else{
                /*
                Si nous ne sommes pas arrivés à l'avant-dernière ligne de la matrice alors nous continuons
                de remplir la matrice avec les coefficients de daubechie. La matrice étant à la base remplie 
                de zéro, nous remplissons chaque ligne avec les 4 coefficients décalés de 2 colonnes à chaque fois.
                Nous rappelons ensuite la fonction afin d'appliquer les mêmes étapes à la ligne suivante de la
                matrice.
                */
                acc2.Change_values(acc1b,acc1,daub_coef_somme);
                acc3.Change_values(acc1b,acc1,daub_coef_soustraction);
                return Daubechie_transform_inverse_recursive(z,acc1+2,acc1b+1,acc2,acc3);
            }
        }
        else{
            // Lorsque acc1b vaut enfin la longueur du signal divisé par 2 alors nous pouvons retourner la matrice de daubechie
            return acc2.ConcatRows2(acc3);
        }
    }

    //Cette fonction permet de calculer la transformée en ondelette inverse
    private Vector Daubechie_transform_inverse(Vector z){
        Vector result;
        int indice_ligne = 0;
        int indice_colonne = 0;

        if(z.VecSize()>2){
            /*
                Si nous sommes dans un cas où la longueur du signal est supérieur à 2 alors nous 
                fixons les paramètres de la fonction daubechie_transform_inverse_recursive et
                nous calculons la matrice de daubechie adéquate au signal. Enfin, nous retournons le
                produit matriciel entre le signal et la matrice de daubechie.
            */
            Matrix A = new(new double[z.VecSize()/2,z.VecSize()]);
            Matrix B = new(new double[z.VecSize()/2,z.VecSize()]);
            Matrix daubechie_matrix = Daubechie_transform_inverse_recursive(z,indice_ligne,indice_colonne,A,B);
            result = z.MatDot(daubechie_matrix);
        }
        else{
            /*
                Si la longueur du signal est égale à 2 nous concaténons le signal avec lui-même puis
                nous répétons les mêmes opérations que dans la boucle précédente sur ce signal concaténé
            */
            Vector z2 = z.VecCon(z).VecCon(z);
            Matrix A = new(new double[z2.VecSize()/2,z2.VecSize()]);
            Matrix B = new(new double[z2.VecSize()/2,z2.VecSize()]);
            Matrix daubechie_matrix = Daubechie_transform_inverse_recursive(z2,indice_ligne,indice_colonne,A,B);
            result = z2.MatDot(daubechie_matrix).SubVec(0,2);
        }

        return result;
    }

    private Vector Smooth_daubechie(Vector signal, int pourcentage_nb){
        // Tout d'abord, nous calculons la transformée en ondelette du signal
        Vector z=Daubechie_transform(signal);
        /*
            Puis nous mettons à 0 tous les éléments du signal transformée dont la valeur absolue 
            est inférieur à un seuil arbitraire
        */
        Vector z_deb=z.MinSignal(0,(int)Math.Ceiling((double)pourcentage_nb*signal.VecSize()/100));
        // Enfin nous retournons le signal débruité après avoir fait une transformée inverse
        return Daubechie_transform_inverse(z_deb);
    }

    // Cette fonction permet de faire une transformée en ondelette multiple. Elle est aussi récursive.
    private Vector Two_level_daubechie_transform(Vector signal,int level,Vector acc){
        // Nous commençons d'abord par faire une première transformée du signal
        Vector z=Daubechie_transform(signal);
        // Puis nous rajoutons le résultat de cette transformée à l'accumulateur
        acc=acc.VecCon(z);
        // Nous regardons si nous avons atteint le niveau maximal
        if(level!=0){
            /*
                Si nous n'avons pas atteint le niveau maximal alors nous rappelons la fonction sur la première moitié
                du signal et ainsi de suite. Cela permet de rajouter à l'accumulateur tous les résultats de ces transformées
            */
            return Two_level_daubechie_transform(z.HalfVec(),level-1,acc);
        }
        else{
            /*
            Lorsque nous avons atteint le niveau maximal alors nous renvoyons acc qui contient les résultats 
            de toutes les transformées successives
            */
            return acc;
        }
    }

    // Cette fonction permet de trier les résultats obtenus avec la fonction Two_level_daubechie_transform
    private Vector[] Multi_daubechie(Vector signal,int level){
        // Tout d'abord, nous appelons la fonction Two_level_daubechie_transform
        Vector z=Two_level_daubechie_transform(signal,level,new([]));
        // Puis nous calculons des puissances de 2
        Vector powers = Vector.PowOf2D(level+1);
        // Nous faisons une somme cumulative de ces puissances de 2
        Vector indexs = powers.CumSum();
        // Enfin nous divisons le tableau z aux indices contenus dans le tableau indexs
        Vector[] z_s = z.Split(indexs);
        // Nous retournons le tableau splitté en plusieurs sous-tableaux
        return z_s;
    }

    /*
        Cette fonction permet de récupérer les éléments du signal transformé plusieurs fois 
        grâce au tableau renvoyé par la fonction multi_daubechie
    */
    private static Vector Multi_trans_daubechie(Vector[] z, int level){
    /*
        Cette fonction permet de sélectionner les éléments
        qui correspondent à la transformée multiple du signal
    */
        z[level-1] = z[level-1].Change_firsts(z[level]);
        if(level != 1){
            return Multi_trans_daubechie(z,level-1);
        }
        else{
            return z[0];
        }
    }

    // Cette fonction permet de faire une transformée multiple à l'aide des fonctions écrites au-dessus
    private Vector Transformee_daubechie_multi(Vector signal, int level) => Multi_trans_daubechie(Multi_daubechie(signal,level),level);

    // Cette fonction permet de splitter le signal transformé plusieurs fois
    private static Vector[] Trans_daubechie_multi_inverse(Vector signal, int level) => signal.Split(Vector.PowOf2A(level+2));

    // Cette fonction permet de faire une transformée multiple inverse et donc de revenir au signal d'origine 
    private Vector Trans_daubechie_inverse(Vector signal, int level, Vector acc, int acc2){
        // Tout d'abord, nous splittons le signal transformé
        Vector[] z=Trans_daubechie_multi_inverse(signal,level);
        // Ensuite, nous rajoutons à acc le bon sous-tableau du signal splitté
        acc=acc.VecCon(z[acc2]);
        // Enfin, nous appliquons une transformée inverse à acc
        acc=Daubechie_transform_inverse(acc);
        if(level != acc2){
            // Si nous n'avons pas atteint le niveau maximal alors nous devons continuer et nous rappelons la fonction
            return Trans_daubechie_inverse(signal,level,acc,acc2+1);
        }
        else{
            // Si nous avons atteint le niveau maximal, alors nous renvoyons le signal d'origine
            return acc;
        }
    }

    // Cette fonction permet de débruiter un signal à l'aide d'une transformée en ondelette multiple
    private Vector Smooth_daubechie_multi(Vector signal,int level,int pourcentage_nb){
        // Tout d'abord, nous calculons la transformée en ondelette multiple du signal
        Vector z=Transformee_daubechie_multi(signal,level);
        /*
            Puis nous mettons à 0 tous les éléments du signal transformée dont la valeur absolue 
            est inférieur à un seuil arbitraire
        */
        Vector z_deb=z.MinSignal(0,(int)Math.Ceiling((double)pourcentage_nb*signal.VecSize()/100));
        // Enfin nous retournons le signal débruité après avoir fait une transformée inverse
        return Trans_daubechie_inverse(z_deb,level,new([]),0);
    }
}