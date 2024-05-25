namespace WaveLibrary;

using ArrayCalc;

public class Haar_Waves // Transformée en ondelette de Haar
{
    // Voici les coefficients de Haar
    private readonly Vector haar_coef_somme = new([1/Math.Sqrt(2),1/Math.Sqrt(2)]);
    private readonly Vector haar_coef_soustraction = new([1/Math.Sqrt(2),-1/Math.Sqrt(2)]);

    private readonly int haar_len;

    public Vector trans_signal;

    public Vector invtrans_signal;

    public Vector smooth_haar;

    public Vector smooth_haar_multi;

    public Haar_Waves(Vector signal, int pct_smooth = 0){
        Console.WriteLine("//////////////////////////////////// Haar ////////////////////////////////////////////");
        haar_len = haar_coef_somme.VecSize();
        trans_signal = W_haar(new([]));
        invtrans_signal = W_haar_inverse(new([]));
        smooth_haar = Smooth_haar(new([]),0);
        smooth_haar_multi = Smooth_haar_multi(new([]),0,0);
        try{
            Console.WriteLine("Début signal");
            signal.VecDisp();
            Console.WriteLine("Fin signal");
            if(signal.VecSize() < 2 || Math.Log2(signal.VecSize()) != (int)Math.Ceiling(Math.Log2(signal.VecSize()))) {
                Console.WriteLine("The size of the vector is not a power of 2 so it does not work");
            }
            else{
                trans_signal = W_haar(signal);
                Console.WriteLine("Début signal transformé");
                trans_signal.VecDisp();
                Console.WriteLine("Fin signal transformé");
                invtrans_signal = W_haar_inverse(trans_signal);
                Console.WriteLine("Début signal transformé inverse");
                invtrans_signal.VecDisp();
                Console.WriteLine("Fin signal transformé inverse");
                smooth_haar = Smooth_haar(signal,pct_smooth);
                Console.WriteLine("Début signal débruité");
                smooth_haar.VecDisp();
                Console.WriteLine("Fin signal débruité");
                if(signal.VecSize() >= 4){
                    smooth_haar_multi = Smooth_haar_multi(signal,(int)Math.Log2(signal.VecSize()-1),pct_smooth);
                    Console.WriteLine("Début signal débruité multi");
                    smooth_haar_multi.VecDisp();
                    Console.WriteLine("Fin signal débruité multi");
                }
                else{
                    Console.WriteLine("The size of the vector is not superior to 4 so the multilevel transformation does not work");
                }
            }
        }
        catch(Exception ex){
            Console.WriteLine(ex.GetType().FullName);
            Console.WriteLine(ex.Message);
            Console.WriteLine(ex.StackTrace);
        }
        Console.WriteLine("//////////////////////////////////////////////////////////////////////////////////////");
    }

    // Cette fonction permet de faire une transformée avec une ondelette de Haar
    private Vector W_haar(Vector signal){
        Vector ret = new([]);
        if(signal.VecSize() > 1 && Math.Log2(signal.VecSize()) == (int)Math.Ceiling(Math.Log2(signal.VecSize()))){
            // Nous divisions le signal en plusieurs pairs de nombres
            Vector[] pairs = signal.Split(Vector.Linspace(signal.VecSize()/haar_len,haar_len));
            // Nous faisons un produit matriciel de ces pairs avec les coefficient de Haar
            Vector tab_somme = new(pairs.Select(v => haar_coef_somme.VecDot(v)).ToArray());
            Vector tab_soustraction = new(pairs.Select(v => haar_coef_soustraction.VecDot(v)).ToArray());
            // Nous renvoyons enfin le résultat de ces produits concaténés
            ret = tab_somme.VecCon(tab_soustraction);
        }
        
        return ret;
    }

    // Cette fonction permet de faire une transformée inverse avec une ondelette de Haar
    private Vector W_haar_inverse(Vector signal){
        Vector ret = new([]);
        if(signal.VecSize() > 1 && Math.Log2(signal.VecSize()) == (int)Math.Ceiling(Math.Log2(signal.VecSize()))){
            // Nous divisons le signal transformée en plusieurs singletons
            Vector[] coefs = signal.Split(Vector.Linspace(signal.VecSize(),1));
            // Nous multiplions chaque singleton de la première moitié du tableau splitté par les coefficients de somme de Haar 
            Vector[] tab_somme_inverse = coefs[0..(signal.VecSize()/2)].Select(v => v[0]*haar_coef_somme).ToArray();
            // Nous multiplions chaque singleton de la seconde moitié du tableau splitté par les coefficients de différence de Haar
            Vector[] tab_soustraction_inverse = coefs[(signal.VecSize()/2)..signal.VecSize()].Select(v => v[0]*haar_coef_soustraction).ToArray();
            // Nous additionnons les deux matrices obtenues
            Vector[] result = Vector.SumTabVec(tab_somme_inverse,tab_soustraction_inverse);
            // Nous pouvons ensuite retourner le signal d'origine
            foreach(Vector v in result){
                ret = ret.VecCon(v);
            }
        }

        return ret;
    }

    private Vector Smooth_haar(Vector signal, int pourcentage_nb){
        Vector ret = new([]);
        if(signal.VecSize() > 1 && pourcentage_nb >= 0 && Math.Log2(signal.VecSize()) == (int)Math.Ceiling(Math.Log2(signal.VecSize()))){
            // Nous faisons d'abord une transformée du signal
            Vector z = W_haar(signal);
            // Puis nous mettons à 0 les valeurs qui sont inférieurs à un certain seuil
            Vector z_deb = z.MinSignal(0,(int)Math.Ceiling((double)pourcentage_nb*signal.VecSize()/100));
            // Enfin nous retournons la transformée inverse du signal débruité
            ret = W_haar_inverse(z_deb);
        }
    
        return ret;
    }

    // Cette fonction permet de faire une transformée en ondelette multiple. Elle est récursive.
    private Vector Two_level_haar_transform(Vector signal,int level, Vector acc){
        // Nous commençons d'abord par faire une première transformée du signal
        Vector z = W_haar(signal);
        // Puis nous rajoutons le résultat de cette transformée à l'accumulateur
        acc = acc.VecCon(z);
        // Nous vérifions si nous avons atteint le niveau maximal
        if(level != 0){
            /*
                Si nous n'avons pas atteint le niveau maximal alors nous rappelons la fonction sur la première moitié
                du signal et ainsi de suite. Cela permet de rajouter à l'accumulateur tous les résultats de ces transformées
            */
            return Two_level_haar_transform(z.HalfVec(),level-1,acc);
        }
        else{
            /*
                Lorsque nous avons atteint le niveau maximal alors nous renvoyons acc qui contient les résultats 
                de toutes les transformées successives
            */
            return acc;
        }
    }

    // Cette fonction permet de trier les résultats obtenus avec la fonction Two_level_haar_transform
    private Vector[] Multi(Vector signal, int level){
        // Tout d'abord, nous appelons la fonction two_level_haar_transform
        Vector z = Two_level_haar_transform(signal, level, new([]));
        // Puis nous calculons des puissances de 2
        Vector powers = Vector.PowOf2D(level + 1);
        // Nous faisons une somme cummulative de ces puissances de 2
        Vector indexs = powers.CumSum();
        // Enfin nous divisons le tableau z aux indices contenus dans le tableau indexs
        Vector[] z_s = z.Split(indexs); 
        // Nous retournons le tableau splitté en plusieurs sous-tableaux
        return z_s;   
    }

    /*
        Cette fonction permet de récupérer les éléments du signal transformé plusieurs fois 
        grâce au tableau renvoyé par la fonction multi
    */
    private static Vector Multi_trans(Vector[] z, int level){
        /*
            Cette fonction permet de sélectionner les éléments
            qui correspondent à la transformée multiple du signal
        */
        z[level-1] = z[level-1].Change_firsts(z[level]);
        if(level != 1){
            return Multi_trans(z,level-1);
        }
        else{
            return z[0];
        }
    }

    // Cette fonction permet de faire une transformée multiple à l'aide des fonctions écrites au-dessus
    private Vector Transformee_multi(Vector signal, int level) => Multi_trans(Multi(signal,level),level);

    // Cette fonction permet de splitter le signal transformé plusieurs fois
    private static Vector[] Trans_multi_inverse(Vector signal, int level) => signal.Split(Vector.PowOf2A(level+2));

    // Cette fonction permet de faire une transformée multiple inverse et donc de revenir au signal d'origine 
    private Vector Trans_inverse(Vector signal, int level, Vector acc, int acc2){
        // Tout d'abord, nous splittons le signal transformé
        Vector[] z=Trans_multi_inverse(signal,level);
        // Ensuite, nous rajoutons à acc le bon sous-tableau du signal splitté
        acc=acc.VecCon(z[acc2]);
        // Enfin, nous appliquons une transformée inverse à acc
        acc=W_haar_inverse(acc);
        if(level != acc2){
            // Si nous n'avons pas atteint le niveau maximal alors nous devons continuer et nous rappelons la fonction
            return Trans_inverse(signal,level,acc,acc2+1);
        }
        else{
            // Si nous avons atteint le niveau maximal, alors nous renvoyons le signal d'origine
            return acc;
        }
    }

    private Vector Smooth_haar_multi(Vector signal,int level,int pourcentage_nb){
        Vector ret = new([]);
        if(signal.VecSize() >= 4 && level > 0 && pourcentage_nb >= 0 && Math.Log2(signal.VecSize()) == (int)Math.Ceiling(Math.Log2(signal.VecSize()))){
            // Tout d'abord, nous calculons la transformée en ondelette multiple du signal
            Vector z=Transformee_multi(signal,level);
            /*
                Puis nous mettons à 0 tous les éléments du signal transformée dont la valeur absolue 
                est inférieur à un seuil arbitraire
            */
            Vector z_deb=z.MinSignal(0,(int)Math.Ceiling((double)pourcentage_nb*signal.VecSize()/100));
            // Enfin nous retournons le signal débruité après avoir fait une transformée inverse
            ret = Trans_inverse(z_deb,level,new([]),0);
        }
        
        return ret;
    }
}