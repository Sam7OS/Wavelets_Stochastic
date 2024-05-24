using WaveLibrary;
using ArrayCalc;

int Min = 100;
int Max = 1000;
Random randNum = new();
double[] tab_d = Enumerable.Repeat(0, 256).Select(i => Min+(Max-Min)*randNum.NextDouble()).ToArray();

Vector vec1 = new(tab_d);

int pct_smooth = 50;

Daub_Waves daub = new(vec1,pct_smooth);

Haar_Waves haar = new(vec1,pct_smooth);

Console.WriteLine("Process done");