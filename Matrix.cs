namespace ArrayCalc;

public class Matrix(double[,] values){

    private readonly double[,] values = values;

    public double this[int row,int column]{
        get { return values[row, column]; }
        set { values[row, column] = value; }
    }

    // Retourne le nombre de lignes
    public int Nb_rows() => values.GetLength(0);

    // Retourne le nombre de colonnes
    public int Nb_cols() => values.GetLength(1);

    // Retourne la colonne numéro num_col
    public Vector GetCol(int num_col) => new(Enumerable.Range(0, Nb_rows()).Select(x => values[x, num_col]).ToArray());

    // Retourne la colonne numéro num_row
    public Vector GetRow(int num_row) => new(Enumerable.Range(0, Nb_cols()).Select(x => values[num_row, x]).ToArray());

    // Permet de modifier les valeurs d'une ligne de la matrice
    public void Change_values(int row, int start, Vector vals){
        if(start >= 0 && start + vals.VecSize() <= Nb_cols()){
            for(int k=start; k<start + vals.VecSize(); k++){
                values[row,k] = vals[k-start];
            }
        }   
    }

    // Retourne la concaténation des lignes de la matrice avec la matrice x ayant le même nombre de lignes
    public Matrix ConcatRows2(Matrix x){
        double[,] concat = new double[2*Nb_rows(),Nb_cols()];
        if(Nb_cols() == x.Nb_cols() && Nb_rows() == x.Nb_rows()){
            for(int k=0;k<Nb_rows();k++){
                for(int h=0;h<Nb_cols();h++){
                    concat[k,h] = this[k,h];
                    concat[Nb_rows()+k,h]=x[k,h];
                }
            }
        }

        return new(concat);
    }

    // Retourne la concaténation des colonnes de la matrice avec la matrice x ayant le même nombre de colonnes
    public Matrix ConcatCols2(Matrix x){
        double[,] concat = new double[Nb_rows(),2*Nb_cols()];
        if(Nb_cols() == x.Nb_cols() && Nb_rows() == x.Nb_rows()){
            for(int k=0;k<Nb_rows();k++){
                for(int h=0;h<Nb_cols();h++){
                    concat[k,h] = this[k,h];
                    concat[k,Nb_cols()+h]=x[k,h];
                }
            }
        }

        return new(concat);
    }    

    // Retourne la concaténation des lignes de la matrice avec la matrice x
    public Matrix ConcatRows(Matrix x){
        double[,] concat = new double[Nb_rows()+x.Nb_rows(),Nb_cols()];
        if(Nb_cols() == x.Nb_cols()){
            for(int k=0;k<Nb_cols();k++){
                for(int h=0;h<Nb_rows();h++){
                    concat[h,k] = this[h,k];
                }
                for(int h=Nb_rows();h<Nb_rows()+x.Nb_rows();h++){
                    concat[h,k] = x[h,k];
                }
            }
        }

        return new(concat);
    }

    // Retourne la concaténation des colonnes de la matrice avec la matrice x
    public Matrix ConcatCols(Matrix x){
        double[,] concat = new double[Nb_rows(),Nb_cols()+x.Nb_cols()];
        if(Nb_rows() == x.Nb_rows()){
            for(int k=0;k<Nb_rows();k++){
                for(int h=0;h<Nb_cols();h++){
                    concat[k,h] = this[k,h];
                }
                for(int h=Nb_cols();h<Nb_cols()+x.Nb_cols();h++){
                    concat[k,h] = x[k,h];
                }
            }
        }

        return new(concat);
    }

    // Permet d'afficher la mtrice sur la console
    public void MatDisp(){
        // Voir pour ToString();
        Console.WriteLine("///////////////////////////////////////////////////////////////////////////////////////////////");
        for(int k=0;k<Nb_rows();k++){
            Console.WriteLine("----------------------------------------------------------------------------------------------");
            Console.Write("| ");
            for(int h=0;h<Nb_cols();h++){
                Console.Write(this[k,h] + " | ");  
            }
            Console.WriteLine();
        }
        Console.WriteLine("///////////////////////////////////////////////////////////////////////////////////////////////");
    }
}