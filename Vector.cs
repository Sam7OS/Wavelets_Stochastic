namespace ArrayCalc;

public class Vector(double[] values)
{
    private double[] values = values;

    public double this[int index]{
        get { return values[index]; }
        set { values[index] = value; }
    }

    // Retourne la taille du vecteur
    public int VecSize() => values.GetLength(0);

    // Retourne le produit scalaire du vecteur avec le vecteur x
    public double VecDot(Vector x) => values.Zip(x.values, (a, b) => a * b).Sum();

    // Retourne la concaténation du vecteur avec le vecteur x
    public Vector VecCon(Vector x) => new([.. values, .. x.values]);

    // Retourne un sous-vecteur entre les indices substart et subend
    public Vector SubVec(int substart, int subend) => new(values[substart..subend]);

    // Retourne le vecteur avec une valeur ajoutée à la fin
    public Vector VecApp(double value) => new([.. values, value]);

    // Retourne le vecteur avec les premières valeurs modifiés
    public Vector Change_firsts(Vector v) => v.VecCon(new(values.Skip(v.VecSize()).ToArray()));

    // Retourne le produit scalaire du vecteur avec une colonne de la matrice x
    public double ColMatDot(Matrix x, int num_col) => VecDot(x.GetCol(num_col));

    // Retourne les produits scalaires du vecteur avec les colonnes de la matrice x
    public Vector MatDot(Matrix x) => new(Enumerable.Range(0,VecSize()).Select(col => ColMatDot(x,col)).ToArray());

    // Retourne la somme cummulative du vecteur 
    public Vector CumSum() => new(values.Select((v,i) => (i == 0)?v:v+values[.. i].Sum()).ToArray());

    // Retourne un vecteur des puissances de 2 croissante
    public static Vector PowOf2A(int length) => new(Enumerable.Range(0,length-1).Select((v,i) => Math.Pow(2,i+1)).ToArray());

    // Retourne un vecteur des puissances de 2 décroissante
    public static Vector PowOf2D(int length) => new(Enumerable.Range(0,length).Select((v,i) => Math.Pow(2,length-i)).ToArray());

    // Retourne un vecteur de taille length partant de 0 avec des multiples de space croissants
    public static Vector Linspace(int length, double space) => new(Enumerable.Range(0,length).Select(v => (v+1)*space).ToArray());

    // Retourne un tableau de sous-vecteurs aux indices contenus dans le vecteur indexs 
    public Vector[] Split(Vector indexs) => indexs.values.Select((v,i) => (i == 0)?SubVec(0,(int)v):SubVec((int)indexs[i-1],(int)v)).ToArray();

    // Pour additionner les valeurs de 2 vecteurs
    public static Vector operator +(Vector x, Vector y) => new(x.values.Zip(y.values, (a, b) => a + b).ToArray());

    // Pour soustraire les valeurs de 2 vecteurs
    public static Vector operator -(Vector x, Vector y) => new(x.values.Zip(y.values, (a, b) => a - b).ToArray());

    // Pour multiplier les valeurs d'un vecteur par un scalaire
    public static Vector operator *(double d, Vector x) => new(x.values.Select(v => v * d).ToArray());

    // Retourne la somme des vecteurs d'un tableau de vecteurs pour chaque indice
    public static Vector[] SumTabVec(Vector[] x, Vector[] y) => x.Zip(y, (a,b) => a + b).ToArray(); 

    private static double ChangeZero(double value,double val_max){
        double ret=value;
        if(ret == val_max + 10){
            ret = 0;
        }
        return ret;
    }

    private static double ChangeMin(double value,double minabs,double val_max){
        double ret = value;
        if(Math.Abs(ret) == minabs){
            ret = val_max + 10;
        }
        return ret;
    }

    private Vector MinRecSignal(int acc,int nb_min,double val_max){
        if(acc == nb_min){
            return new(values.Select(v => ChangeZero(v,val_max)).ToArray());
        }
        else{
            double minabs = values.Select(Math.Abs).Min();
            values = values.Select(v => ChangeMin(v,minabs,val_max)).ToArray();
            return MinRecSignal(acc+1,nb_min,val_max);
        }
    }

    // Cette fonction permet de mettre un nombre renseigné de minimums à 0 
    public Vector MinSignal(int acc,int nb_min){
        double val_max=values.Select(Math.Abs).Max();
        int nb_duplic = VecSize() - values.Select(Math.Abs).Distinct().Count();
        return MinRecSignal(acc,nb_min+nb_duplic,val_max);
    }

    // Retourne la première moitié du vecteur
    public Vector HalfVec(){
        return new(values.Chunk(VecSize()/2).ToArray()[0]);
    }

    // Affiche le vecteur sur la console
    public void VecDisp(){
        Console.Write("| ");
        Array.ForEach(values, v => {Console.Write(v + " | ");});
        Console.WriteLine();
    }

    // Implémentation de la méthode ToString
    public override string ToString()
    {
        string str = "| ";
        _ = values.All(v => { str += v.ToString() + " | "; return true; });
        return str;
    }
}