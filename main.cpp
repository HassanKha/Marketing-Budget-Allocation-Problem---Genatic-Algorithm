#include <bits/stdc++.h>
using namespace std;

class Fittnesscalculation
{

private:
    vector<double>ROI;


public:
    Fittnesscalculation() {}
    Fittnesscalculation(vector<double>ROII)
    {
        ROI=ROII;
    }

    double getGenesFittness(vector<double>genes)
    {
        //  string gene=c.getGenes();
        double fittnes=0;
        for(int i=0; i<genes.size(); ++i)
        {
            double pro=(genes[i]*ROI[i])/100.0;
            fittnes+=pro;
        }
        return fittnes;
    }

};


class Chromosome
{

private:
    vector<double>genes;
    vector<pair<double,double>>constains;
    double genesFitness;

public:

    Chromosome()
    {
    }

    Chromosome(vector<pair<double,double>>constainss)
    {

        constains=constainss;
    }
    void generateGenesRandom(double totalBudget,int chaneelSizes)
    {
        double totalLower=getTotalLower();
        double currBudget;
        double upperConstrain;
        for(int i=0; i<chaneelSizes; ++i)
        {
            upperConstrain=constains[i].second;
            totalLower=abs(totalLower-constains[i].first);
            currBudget=abs(totalBudget-totalLower);
            if(constains[i].second>currBudget)
                upperConstrain=currBudget;
            double randomBudgetInvestInChannel=fRand(constains[i].first,upperConstrain);
           // cout<< "what is it :" << randomBudgetInvestInChannel <<endl;
            genes.push_back(randomBudgetInvestInChannel);
            double reamingOCurrBudget=currBudget-randomBudgetInvestInChannel;
            totalBudget=reamingOCurrBudget+totalLower;


        }


   //PrintGenes();


    }
    double GetSumofGenes(){
        double Sum=0;
        for(int i=0 ; i<genes.size() ; i++){
            Sum+=genes[i];
        }
        return Sum;
    }

void setConstrains(  vector<pair<double,double>>constainss){
constains=constainss;
}
    void PrintGenes(){
        double k=0;
    for(int i=0 ; i<genes.size() ; i++){
        cout << i+1 << "th Channel ---->>> "<< genes[i] << endl ;
        k+=genes[i];
    }
  //  cout <<endl;
cout << "Cost : " << k <<endl;
cout <<endl;
    }
    int getGenesSize()
    {
        return genes.size();
    }
    double getGenesFitt(vector<double>ROI,vector<double>g)
    {
        Fittnesscalculation FC(ROI);
        genesFitness=FC.getGenesFittness(g);
        return genesFitness;

    }
    double getLowerBound(int i)
    {
        return constains[i].first;

    }
    double getUpperBound(int i)
    {
        return constains[i].second;

    }
    double getGene(int ind)
    {
        return genes[ind];

    }
    void setGene(double c,int ind)
    {
        genes[ind]=c;

    }
    void PushRange(double lowerbound,double Upperbound)
    {

        constains.push_back({lowerbound,Upperbound});

    }
    vector<double> getGenes()
    {
        return genes;
    }

    void setChromosome(Chromosome C)
    {

        genes=C.getGenes();
    }
    void setGenes(vector<double>G)
    {
        genes=G;
    }
    void PushGene(double Gene){
    genes.push_back(Gene);
    }

    double getTotalLower()
    {
        double sum=0;
        for(int i=0; i<constains.size(); ++i)
        {
            sum+=constains[i].first;
        }

        return sum;

    }

    double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }

};




class Generations
{
private:
    vector<Chromosome>chromosome;
    int  current_generation;
    int MaxofGeneration;
vector<double>ROI;


public:
    Generations(int n,bool check,double totalBudget,int chaneelSizes,vector<pair<double,double>>Range)
    {

        if(check)
        {
            current_generation=1;
            creatRandomGeneration(n,totalBudget,chaneelSizes,Range);
        }
    }
    void creatRandomGeneration(int n,double totalBudget,int chaneelSizes,vector<pair<double,double>>Range)
    {
        for(int i=0; i<n; ++i)
        {
            Chromosome NC(Range);

            NC.generateGenesRandom(totalBudget,chaneelSizes);
            // NC.PrintGenes();
            chromosome.push_back(NC);

        }
    }
    void PrintChromosomes(){
        cout <<endl;
    for(int i=0; i <chromosome.size();i++){

          chromosome[i].PrintGenes();
    }
cout <<endl;
    }
    void SetMaxGenerations(int MG)
    {
        MaxofGeneration=MG;
    }
    int GetMaxGenerations()
    {
        return MaxofGeneration;
    }
    int bestChromosomeindex(vector<double>ROI)
    {
        Chromosome chromosomeFitt;
        chromosomeFitt.setChromosome(chromosome[0]);
int in=0;
        for(int i=0; i<chromosome.size(); ++i)
        {

            if(chromosomeFitt.getGenesFitt(ROI,chromosomeFitt.getGenes())<chromosome[i].getGenesFitt(ROI,chromosome[i].getGenes()))
            {
                chromosomeFitt.setChromosome(chromosome[i]);
                in=i;
            }

        }
        return in;

    }
    Chromosome bestChromosome(vector<double>ROI)
    {
        Chromosome chromosomeFitt;
        chromosomeFitt.setChromosome(chromosome[0]);

        for(int i=0; i<chromosome.size(); ++i)
        {

            if(chromosomeFitt.getGenesFitt(ROI,chromosomeFitt.getGenes())<chromosome[i].getGenesFitt(ROI,chromosome[i].getGenes()))
            {
                chromosomeFitt.setChromosome(chromosome[i]);
            }

        }
        return chromosomeFitt;

    }

    int worstChromosome(vector<double>ROI)
    {
        Chromosome chromosomeFitt;
        chromosomeFitt.setChromosome(chromosome[0]);

      int index=0;
        for(int i=0; i<chromosome.size(); ++i)
        {
            if(chromosomeFitt.getGenesFitt(ROI,chromosomeFitt.getGenes())>chromosome[i].getGenesFitt(ROI,chromosome[i].getGenes()))
            {
                chromosomeFitt.setChromosome(chromosome[i]);
                index=i;
            }

        }
       // cout << "worst" <<endl;
//cout << chromosomeFitt.getGenesFitt(ROI,chromosomeFitt.getGenes());
//cout << "worst" <<endl;
        return index;

    }
    double generationFitness(vector<double>ROI)
    {
        double fitnesssum=0;
        for(int i=0; i<chromosome.size(); ++i)
        {
            fitnesssum+=chromosome[i].getGenesFitt(ROI,chromosome[i].getGenes());
        }
        return fitnesssum;
    }

    int GenerationSize()
    {
        return chromosome.size();
    }
     void SetChromosomes(vector<Chromosome>V)
    {
       chromosome=V;
    }

    Chromosome getChromosomeIndex(int n)
    {
        return chromosome[n];
    }
void setChromo(Chromosome c,int i){
chromosome[i]=c;
}
    vector<Chromosome>getChromosome()
    {
        return chromosome;
    }

    void pushChromosome(Chromosome C)
    {
        chromosome.push_back(C);
    }

    void upgrade_genration(int LastGeneration)
    {
        current_generation=LastGeneration+1;
    }

    int getCurrentGeneration()
    {
        return current_generation;
    }
   void RemoveChromo(int k){
       vector<Chromosome>CC;
       int x=0;
   for(int i=0 ; i<chromosome.size(); i++){
    if(i!=k){
        CC.push_back(chromosome[i]);

    }
   }
   SetChromosomes(CC);

   }


};



class GA
{
private:

    double mutationRate = 0.5;
    int tournamentSize = 2;


public:
    Generations evolvePopulation(Generations generation,vector<double>ROI,double totalBudget,int chaneelSizes,vector<pair<double,double>>Range)
    {
        Generations newGeneration(generation.GenerationSize(), false,totalBudget,chaneelSizes,Range);

        newGeneration.upgrade_genration(generation.getCurrentGeneration());


    Generations DG=generation;

        for (int i = 0; i < generation.GenerationSize(); i++)
        {
             pair<Chromosome,int> c1 = tournamentSelection(DG,ROI,generation,totalBudget, chaneelSizes,Range);



          c1.first.setConstrains(Range);

            newGeneration.pushChromosome(c1.first);
            DG.RemoveChromo(c1.second);


        }


for (int i = 0; i < newGeneration.GenerationSize(); i+=2)
        {



             vector<Chromosome> newc = crossover(newGeneration.getChromosomeIndex(i), newGeneration.getChromosomeIndex(i+1),totalBudget);

        newc[0].setConstrains(Range);
           newc[1].setConstrains(Range);


              newGeneration.setChromo(newc[0],i);
                newGeneration.setChromo(newc[1],i+1);
        }


        for (int i = 0; i < newGeneration.GenerationSize(); i++)
        {


    Chromosome m= NoNUniformMutate(newGeneration.getChromosomeIndex(i),newGeneration,totalBudget);

        //   Chromosome m= UniformMutate(newGeneration.getChromosomeIndex(i),totalBudget);

          m.setConstrains(Range);
          newGeneration.setChromo(m,i);


        }

      newGeneration= Replacement(generation,newGeneration,ROI,Range);


        return newGeneration;
    }

Generations Replacement(Generations generation,Generations newGeneration,vector<double>ROI,vector<pair<double,double>>Range){
int i= newGeneration.worstChromosome(ROI);
Chromosome F= generation.bestChromosome(ROI);

F.setConstrains(Range);

newGeneration.setChromo(F,i);
Chromosome F1= newGeneration.bestChromosome(ROI);

return newGeneration;
}

    pair<Chromosome,int> tournamentSelection(Generations generation,vector<double>ROI,Generations generationn,double totalBudget,int chaneelSizes,vector<pair<double,double>>Range)
    {

        Generations Generationtournament(tournamentSize, false,totalBudget, chaneelSizes,Range);


        for (int i = 0; i < tournamentSize; i++)
        {

         int R=rand()%generation.GenerationSize();

            Generationtournament.pushChromosome(generation.getChromosomeIndex(R));

        }

        Chromosome fittest = Generationtournament.bestChromosome(ROI);
        int h=generation.bestChromosomeindex(ROI);
pair<Chromosome,int>b;
b.first=fittest;
b.second=h;
        return b;
    }
 vector<Chromosome> crossover(Chromosome c1, Chromosome c2,double totalBudget)
    {
vector<Chromosome>Children;
vector<pair<int,int>>Points;
    for(int i=0 ; i<c1.getGenes().size() ; i++){
        for(int j=i+1 ; j<c2.getGenes().size() ; j++){
        Points.push_back({i,j});

    }

    }

    vector<pair<int,int>>AvPoints;
    for(int i=0 ; i<Points.size();i++){

        int s=c1.getGenes().size();
        int ind1=Points[i].first;
        int ind2=Points[i].second;

        vector<double>v=c1.getGenes();
        vector<double>v2=c2.getGenes();

        if(ind1>ind2)
        {

            for(int i=ind2; i<ind1; ++i)
            {

                double tempc1=v[i];
                v[i]=v2[i];
                v2[i]=tempc1;

            }
        }
        else
        {

            for(int i=ind1; i<ind2; ++i)
            {
                double tempc1=v[i];
                v[i]=v2[i];
                v2[i]=tempc1;
            }

        }

        vector<Chromosome>ans;
       Chromosome nc,nc2;
       nc.setGenes(v);
       nc2.setGenes(v2);

       ans.push_back(nc);
       ans.push_back(nc2);

if(nc.GetSumofGenes()<=totalBudget && nc2.GetSumofGenes()<=totalBudget){
    AvPoints.push_back({ind1,ind2});
//cout << ind1 << " " << ind2 <<endl;
}

    }

    if(AvPoints.size()==0){
        Children.push_back(c1);
        Children.push_back(c2);
       // cout <<"no" <<endl;
        return Children;
    }
    else {
            // cout <<"yes" <<endl;
        int ind=rand()%AvPoints.size();
        int ind1 = AvPoints[ind].first;
        int ind2 = AvPoints[ind].second;

        vector<double>v=c1.getGenes();
        vector<double>v2=c2.getGenes();

        if(ind1>ind2)
        {

            for(int i=ind2; i<ind1; ++i)
            {

                double tempc1=v[i];
                v[i]=v2[i];
                v2[i]=tempc1;

            }
        }
        else
        {

            for(int i=ind1; i<ind2; ++i)
            {
                double tempc1=v[i];
                v[i]=v2[i];
                v2[i]=tempc1;
            }

        }

        vector<Chromosome>ans;
       Chromosome nc,nc2;
       nc.setGenes(v);
       nc2.setGenes(v2);

       ans.push_back(nc);
       ans.push_back(nc2);






       return ans;

    }


    }

    double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }

    Chromosome UniformMutate(Chromosome c,double cost)
    {
vector<double>Genes;

        for (int i = 0; i < c.getGenes().size(); i++)
        {
            double Random1 = ((double) rand() / (RAND_MAX)) ;
            vector<double>Mychromosome=c.getGenes();
            double Delta;

            if (Random1 <= mutationRate)
            {

                double DeltaLower=Mychromosome[i]-c.getLowerBound(i);
                Delta = DeltaLower;

                double Random2=fRand(0,Delta);

                double newGene = Mychromosome[i] - Random2;

Genes.push_back(newGene);
            }
            else
            {

                double DeltaUpper = c.getUpperBound(i)-Mychromosome[i];
                Delta = DeltaUpper;
                double Random2 = fRand(0,Delta);
                double newGene = Mychromosome[i] + Random2;

Genes.push_back(newGene);
            }
        }

        c.setGenes(Genes);

while(c.GetSumofGenes()>cost){
    vector<double>Genes;

        for (int i = 0; i < c.getGenes().size(); i++)
        {
            double Random1 = ((double) rand() / (RAND_MAX)) ;
            vector<double>Mychromosome=c.getGenes();
            double Delta;

            if (Random1 <= mutationRate)
            {

                double DeltaLower=Mychromosome[i]-c.getLowerBound(i);
                Delta = DeltaLower;

                double Random2=fRand(0,Delta);

                double newGene = Mychromosome[i] - Random2;

Genes.push_back(newGene);
            }
            else
            {

                double DeltaUpper = c.getUpperBound(i)-Mychromosome[i];
                Delta = DeltaUpper;
                double Random2 = fRand(0,Delta);
                double newGene = Mychromosome[i] + Random2;

Genes.push_back(newGene);
            }
        }
        c.setGenes(Genes);
}


return c;

    }
    Chromosome NoNUniformMutate(Chromosome c,Generations generation,double cost)
    {

        for (int i = 0; i < c.getGenesSize(); i++)
        {
            double Random1 = ((double) rand() / (RAND_MAX)) ;
            vector<double>Mychromosome=c.getGenes();
            double Delta;
            if (Random1 <= mutationRate)
            {
                double DeltaLower=Mychromosome[i]-c.getLowerBound(i);
                Delta = DeltaLower;

                double r=rand()%2;
                double b=rand()%6;
                double rule = pow((1-(generation.getCurrentGeneration()/generation.GetMaxGenerations())),b);
                double NewDelta= Delta * (1-pow(r,rule));
                double newGene = Mychromosome[i] - NewDelta;
                c.setGene(newGene,i);
            }
            else
            {
                double DeltaUpper = c.getUpperBound(i)-Mychromosome[i];
                Delta = DeltaUpper;

                double r=rand()%2;
                double b=rand()%6;
                double rule = pow((1-(generation.getCurrentGeneration()/generation.GetMaxGenerations())),b);
                double NewDelta= Delta * (1-pow(r,rule));
                double newGene = Mychromosome[i] + NewDelta;
                c.setGene(newGene,i);
            }
        }

        while(c.GetSumofGenes()>cost){
            for (int i = 0; i < c.getGenesSize(); i++)
        {
            double Random1 = ((double) rand() / (RAND_MAX)) ;
            vector<double>Mychromosome=c.getGenes();
            double Delta;
            if (Random1 <= mutationRate)
            {
                double DeltaLower=Mychromosome[i]-c.getLowerBound(i);
                Delta = DeltaLower;

                double r=rand()%2;
                double b=rand()%6;
                double rule = pow((1-(generation.getCurrentGeneration()/generation.GetMaxGenerations())),b);
                double NewDelta= Delta * (1-pow(r,rule));
                double newGene = Mychromosome[i] - NewDelta;
                c.setGene(newGene,i);
            }
            else
            {
                double DeltaUpper = c.getUpperBound(i)-Mychromosome[i];
                Delta = DeltaUpper;

                double r=rand()%2;
                double b=rand()%6;
                double rule = pow((1-(generation.getCurrentGeneration()/generation.GetMaxGenerations())),b);
                double NewDelta= Delta * (1-pow(r,rule));
                double newGene = Mychromosome[i] + NewDelta;
                c.setGene(newGene,i);
            }
        }

        }
return c ;
    }

};


int main()
{
    srand(time(0));

cout << " Enter the marketing budget (in thousands): " <<endl;
int Cost;
cin >> Cost ;
cout << " Enter the number of marketing channels: " <<endl;
int Channels;
cin >> Channels ;
cout << "Enter ROI (in %) of each channel separated by space:" <<endl;
vector<double>ROI;

int roi;
for(int i=0 ; i<Channels ; i++){
    cin >> roi ;
    ROI.push_back(roi);
}
vector<pair<double,double>>Range;
string lower,upper;
cout << "Enter the lower (k) and upper bounds (%) of investment in each channel:(enter x if there is no bound)"<<endl;
for(int i=0 ; i<Channels ; i++){
    cin >> lower >> upper ;
    if(lower=="x"){
             std::string s = std::to_string(0);
        lower=s;
    }
     if(upper=="x"){
            std::string s = std::to_string(Cost);
        lower=s;
    }
       double Lower = std::stod(lower);
       double Upper = std::stod(upper);
       Upper=(Cost*Upper)/100;
       Range.push_back({Lower,Upper});
}
int loop=20;
ofstream myfile;
 myfile.open ("NonMutate.txt");
 Generations MyAnswers(loop,false,Cost,Channels,Range);
while(loop--){
 Generations G (10,true,Cost,Channels,Range);
    GA ga ;
 int GenerattionsSize=10;
Generations BestGenerations=G;



while(GenerattionsSize--)
    {

        Generations NewGeneration = ga.evolvePopulation(BestGenerations,ROI,Cost,Channels,Range);


        if(NewGeneration.generationFitness(ROI)>= BestGenerations.generationFitness(ROI))
        {
            BestGenerations = NewGeneration;
        }
        else
        {
            break;
        }


    }
double Answer = BestGenerations.bestChromosome(ROI).getGenesFitt(ROI,BestGenerations.bestChromosome(ROI).getGenes());
cout << "The total profit : " << BestGenerations.bestChromosome(ROI).getGenesFitt(ROI,BestGenerations.bestChromosome(ROI).getGenes()) <<endl;
cout << "The Best Chromosome : " <<endl;
BestGenerations.bestChromosome(ROI).PrintGenes();
vector <double> Chromo = BestGenerations.bestChromosome(ROI).getGenes();
MyAnswers.pushChromosome(BestGenerations.bestChromosome(ROI));

 myfile << endl;
myfile << "The total profit : ";
  myfile << Answer;
  myfile << endl;
for(int i=0 ; i<Chromo.size() ; i++){
         myfile << i+1;
myfile << "th Channel ----->>> ";
    myfile << Chromo[i];
     myfile << endl;
}
 myfile << endl;
}
myfile << "------------------------------------------------------------------";
myfile << endl;
double BestAnswer = MyAnswers.bestChromosome(ROI).getGenesFitt(ROI,MyAnswers.bestChromosome(ROI).getGenes());
myfile << "Best Total Profit : ";
myfile << BestAnswer;
myfile << endl;
vector <double> Chromoo = MyAnswers.bestChromosome(ROI).getGenes();
for(int i=0 ; i<Chromoo.size() ; i++){
         myfile << i+1;
myfile << "th Channel ----->>> ";
    myfile << Chromoo[i];
     myfile << endl;
}
 myfile << endl;


 myfile.close();

    return 0;
}
/*
8
12
7
11
*/

/*
2.7 58
20.5 100
0 18
10 100
*/
