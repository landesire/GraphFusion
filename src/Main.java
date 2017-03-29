import java.io.*;
import java.util.*;


public class Main {
    static  int GraphNumber;
    static  int PopulationSize; //种群数量最好为偶数
    static  int generationNumber;
    static  int ElitismNumber;
    static  int DiversityThreshold;
    public static void main(String[] args) throws Exception{
        GraphNumber=Integer.parseInt(args[0]);
        PopulationSize=Integer.parseInt(args[1]);
        generationNumber=Integer.parseInt(args[2]);
        ElitismNumber = PopulationSize/50;
        DiversityThreshold=PopulationSize/50;
        String GraphDatafilename=args[3];
        String fitnessOuptutfilename=args[4];
        Writer w=new FileWriter(fitnessOuptutfilename,true);
        for(int count=0;count<100;count++)
        {
            GeneticAlgo(GraphDatafilename,fitnessOuptutfilename,w);
        }

        w.close();

    }


    public static void  GeneticAlgo (String GraphDatafilename,String fitnessOuptutfilename, Writer w)throws Exception
    {
        int Graphset[][][]=new int[GraphNumber][][];

        int fitness[];
        int VertexNumberAfterFusion[]=new int[PopulationSize];
        int EdgeNumberAfterFusion[]=new int[PopulationSize];

        boolean VisitedNumberSet[][];


        int firstPopulationMatrix[][][];
        int index,col,i,j,populationIndex;
        int bestIndividualIndex=-1;


        List<List<List<Integer>>> firstPopulation;
        List<List<List<Integer>>> nextGenerationParent;
        List<List<List<Integer>>> nextGenerationChild;
        Queue<Individual> ElitismPQ;
        long min_fitness=Long.MAX_VALUE;


        //计算所有图的所有节点数
        int VertexALLNumber=0;

        //计算所有图的之前节点数
        int VertexPreNumber[]=new int[GraphNumber];

        //保持种群多样性时的访问数组
        boolean DiversityVisitedPopulationSet[]=new boolean[PopulationSize];
        Arrays.fill(DiversityVisitedPopulationSet,false);




        //生成初始图集合
        GenerateGraph(GraphDatafilename,Graphset);

        for(index=0;index<GraphNumber;index++)
        {
            VertexPreNumber[index]=VertexALLNumber;
            // System.out.println("VertextPreNumber"+index+": "+VertexPreNumber[index]);
            VertexALLNumber=VertexALLNumber+Graphset[index].length;
        }
        //System.out.println(VertexALLNumber);
        // System.out.println(VertexPreNumber[0]);
        //System.out.println(VertexPreNumber[1]);





        //生成初始图矩阵
        firstPopulationMatrix=generateFirstPopulaitonMatrix(Graphset,VertexALLNumber);
        //System.out.println(Arrays.deepToString(firstPopulationMatrix[0]));
        // System.out.println(Arrays.deepToString(firstPopulationMatrix[1]));

        //生成初始下三角矩阵编码
        firstPopulation=FullMatrixToLadderMatrix(firstPopulationMatrix,VertexPreNumber,Graphset,VertexALLNumber);

        //测试初始下三角矩阵编码
        /*
        System.out.println("测试初始下三角矩阵编码~~~~~~~~~~~~~~~~~~");
        System.out.println(firstPopulation.get(0));
        System.out.println(firstPopulation.get(1));
        System.out.println(firstPopulation.get(2));

        System.out.println("测试初始下三角矩阵编码结束~~~~~~~~~~~~~~~~~~");
        */

        //计算适应度和冲突数组
        VisitedNumberSet=caculateVistitedNumberSet(firstPopulation,VertexALLNumber,VertexPreNumber);
        //测试VisitedNumberSet的计算
        //System.out.println("测试VisitedNumberSet的计算是否正确~~~~~~~~~~~~~~~~~~");
        //System.out.println(Arrays.toString(VisitedNumberSet[0]));
        //System.out.println(Arrays.toString(VisitedNumberSet[1]));
        //System.out.println(Arrays.toString(VisitedNumberSet[2]));
        //System.out.println("测试VisitedNumberSet的计算是否正确结束~~~~~~~~~~~~~~~~~~");

        VertexNumberAfterFusion=CaculateVertexNumberAfterFusion(firstPopulation,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        EdgeNumberAfterFusion=CaculateEdgeNumberAfterFusion(firstPopulation,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        fitness=CaculateFitness(firstPopulation,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        //min_fitness=Long.MAX_VALUE;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            if(fitness[populationIndex]<min_fitness)
            {
                min_fitness=fitness[populationIndex];
                bestIndividualIndex=populationIndex;
            }
        }
        System.out.println("对于初始化种群的最优个体的index是: "+bestIndividualIndex+" "+min_fitness);
        //w.write("初始最优个体点数 "+VertexNumberAfterFusion[bestIndividualIndex]+"\n");
        //w.write("初始最优个体边数 "+EdgeNumberAfterFusion[bestIndividualIndex]+"\n");
        //w.write("初始最优个体"+" "+min_fitness+"\n");
        //初代最优个体fitness

        w.write(min_fitness+" ");
        //System.out.println("测试fitness的计算是否正确,包含点和边~~~~~~~~~~~~~~~~~~~~");
        //System.out.println(Arrays.toString(VertexNumberAfterFusion));
        //System.out.println(Arrays.toString(EdgeNumberAfterFusion));
        //System.out.println(Arrays.toString(fitness));
        //System.out.println("测试fitness的计算是否正确,包含点和边,测试结束~~~~~~~~~~~~");

        //精英主义,保留fitness最小的N个个体
        ElitismPQ=Elitism(firstPopulation,fitness,VertexALLNumber,VisitedNumberSet);

        //按照适应度倒数轮盘选择
        nextGenerationParent=Select(fitness,firstPopulation);

        //选好Parent之后,VisitedNumberSet需要对应的更新
        VisitedNumberSet=caculateVistitedNumberSet(nextGenerationParent,VertexALLNumber,VertexPreNumber);
        /*
        System.out.println("测试经过select操作以后的编码是否合理~~~~~~~~~~~~~~~~~~~~");
        System.out.println(nextGenerationParent.get(0));
        System.out.println(nextGenerationParent.get(1));
        System.out.println(nextGenerationParent.get(2));

        System.out.println(Arrays.toString(VisitedNumberSet[0]));
        System.out.println(Arrays.toString(VisitedNumberSet[1]));
        System.out.println(Arrays.toString(VisitedNumberSet[2]));
        System.out.println("测试经过select操作以后的编码是否合理结束~~~~~~~~~~~~~~~~~~~~");
        */
        //System.out.println("crossover befor VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
        //System.out.println("crossover befor VisitedSet"+Arrays.toString(VisitedNumberSet[1]));
        //粗粒度交叉
        /*
        System.out.println("测试经过crossover操作编码是否合理~~~~~~~~~~~~~~~~~~~~");
        System.out.println(nextGenerationParent.get(0));
        System.out.println(nextGenerationParent.get(1));
        System.out.println(nextGenerationParent.get(2));
        System.out.println(nextGenerationParent.get(3));
        System.out.println(nextGenerationParent.get(4));
        System.out.println(nextGenerationParent.get(5));
        System.out.println(Arrays.toString(VisitedNumberSet[0]));
        System.out.println(Arrays.toString(VisitedNumberSet[1]));
        System.out.println(Arrays.toString(VisitedNumberSet[2]));
        System.out.println(Arrays.toString(VisitedNumberSet[3]));
        System.out.println(Arrays.toString(VisitedNumberSet[4]));
        System.out.println(Arrays.toString(VisitedNumberSet[5]));
        */
        nextGenerationChild=coarse_grained_crossover(nextGenerationParent,VertexPreNumber,VisitedNumberSet,VertexALLNumber);
        /*
        System.out.println(nextGenerationChild.get(0));
        System.out.println(nextGenerationChild.get(1));
        System.out.println(nextGenerationChild.get(2));
        System.out.println(nextGenerationChild.get(3));
        System.out.println(nextGenerationChild.get(4));
        System.out.println(nextGenerationChild.get(5));
        System.out.println(Arrays.toString(VisitedNumberSet[0]));
        System.out.println(Arrays.toString(VisitedNumberSet[1]));
        System.out.println(Arrays.toString(VisitedNumberSet[2]));
        System.out.println(Arrays.toString(VisitedNumberSet[3]));
        System.out.println(Arrays.toString(VisitedNumberSet[4]));
        System.out.println(Arrays.toString(VisitedNumberSet[5]));
        System.out.println("测试经过crossover操作编码是否合理结束~~~~~~~~~~~~~~~~~~~~");
        */

        //交叉后对全体群体重新计算VisitedNumberSet
        //VisitedNumberSet=caculateVistitedNumberSet(nextGenerationChild,VertexALLNumber,VertexPreNumber);
        //System.out.println(nextGenerationChild.get(0));
        //System.out.println(nextGenerationChild.get(1));

        //System.out.println("crossover after VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
        //System.out.println("crossover after VisitedSet"+Arrays.toString(VisitedNumberSet[1]));

        //细粒度变异
        //System.out.println("mutation befor VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
        //System.out.println("mutation befor VisitedSet"+Arrays.toString(VisitedNumberSet[1]));
        /*
        System.out.println("测试经过mutation操作编码是否合理~~~~~~~~~~~~~~~~~~~~");
        System.out.println(nextGenerationChild.get(0));
        System.out.println(nextGenerationChild.get(1));
        System.out.println(nextGenerationChild.get(2));
        System.out.println(Arrays.toString(VisitedNumberSet[0]));
        System.out.println(Arrays.toString(VisitedNumberSet[1]));
        System.out.println(Arrays.toString(VisitedNumberSet[2]));
        */
        fine_grained_mutation(nextGenerationChild,Graphset,VertexALLNumber,VertexPreNumber,VisitedNumberSet);
        /*
        System.out.println(nextGenerationChild.get(0));
        System.out.println(nextGenerationChild.get(1));
        System.out.println(nextGenerationChild.get(2));
        System.out.println(Arrays.toString(VisitedNumberSet[0]));
        System.out.println(Arrays.toString(VisitedNumberSet[1]));
        System.out.println(Arrays.toString(VisitedNumberSet[2]));
        System.out.println("测试经过mutation操作编码是否合理结束~~~~~~~~~~~~~~~~~~~~");
        */
        //变异操作对于visitedNumberSet的操作已经包含在内,不需要重新计算

        VisitedNumberSet=caculateVistitedNumberSet(nextGenerationChild,VertexALLNumber,VertexPreNumber);
        fitness=CaculateFitness(nextGenerationChild,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        //System.out.println("mutation after VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
        //System.out.println("mutation after VisitedSet"+Arrays.toString(VisitedNumberSet[1]));
        //利用精英终于生成的E个个体和原来的N个个体混合,取前N个最优的
        // System.out.println("精英主义替换开始~~~~~~~~~~~~~~~~~~~");
        // System.out.println(nextGenerationChild);
        //System.out.println("精英主义替换前的fitness~");
        ElitismReplace(nextGenerationChild,ElitismPQ,fitness,VisitedNumberSet,VertexALLNumber);
        //System.out.println(nextGenerationChild);
        //System.out.println("精英主义替换完成~~~~~~~~~~~~~~~~~~~");

        keepDiversty(nextGenerationChild,DiversityVisitedPopulationSet,Graphset,VertexALLNumber,VertexPreNumber);
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            if(fitness[populationIndex]<min_fitness)
            {
                min_fitness=fitness[populationIndex];
            }
        }
        //w.write(0+" "+min_fitness+"\n");
        for(int count=1;count<generationNumber;count++)
        {
            //每轮补充新个体,变异个体以后要重新计算VisitedNumberSet(冲突数组)
            VisitedNumberSet=caculateVistitedNumberSet(nextGenerationChild,VertexALLNumber,VertexPreNumber);
            fitness=CaculateFitness(nextGenerationChild,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
            //System.out.println("fitness arrays~"+Arrays.toString(fitness));
            ElitismPQ=Elitism(nextGenerationChild,fitness,VertexALLNumber,VisitedNumberSet);
            for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
            {
                if(fitness[populationIndex]<min_fitness)
                {
                    min_fitness=fitness[populationIndex];
                }
            }
            nextGenerationParent=Select(fitness,nextGenerationChild);
            //选好Parent之后,VisitedNumberSet需要对应的更新
            VisitedNumberSet=caculateVistitedNumberSet(nextGenerationParent,VertexALLNumber,VertexPreNumber);

            //System.out.println(nextGenerationParent.get(0));
            //System.out.println(nextGenerationParent.get(1));

            // System.out.println("crossover befor VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
            //System.out.println("crossover befor VisitedSet"+Arrays.toString(VisitedNumberSet[1]));
            nextGenerationChild=coarse_grained_crossover(nextGenerationParent,VertexPreNumber,VisitedNumberSet,VertexALLNumber);

            //VisitedNumberSet=caculateVistitedNumberSet(nextGenerationChild,VertexALLNumber,VertexPreNumber);
            //System.out.println(nextGenerationChild.get(0));
            //System.out.println(nextGenerationChild.get(1));

            //System.out.println("crossover after VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
            //System.out.println("crossover after VisitedSet"+Arrays.toString(VisitedNumberSet[1]));



            fine_grained_mutation(nextGenerationChild,Graphset,VertexALLNumber,VertexPreNumber,VisitedNumberSet);
            //System.out.println("精英主义替换开始~~~~~~~~~~~~~~~~~~~");
            //System.out.println(nextGenerationChild);
            //System.out.println("精英主义替换前的fitness~");
            VisitedNumberSet=caculateVistitedNumberSet(nextGenerationChild,VertexALLNumber,VertexPreNumber);
            fitness=CaculateFitness(nextGenerationChild,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);

            ElitismReplace(nextGenerationChild,ElitismPQ,fitness,VisitedNumberSet,VertexALLNumber);
            //System.out.println(nextGenerationChild);
            //System.out.println("精英主义替换完成~~~~~~~~~~~~~~~~~~~");

            Arrays.fill(DiversityVisitedPopulationSet,false);
            keepDiversty(nextGenerationChild,DiversityVisitedPopulationSet,Graphset,VertexALLNumber,VertexPreNumber);

            // nextGenerationChild.get(PopulationSize-1).clear();
            //nextGenerationChild.get(PopulationSize-2).clear();
            //nextGenerationChild.get(PopulationSize-1).addAll(ElitismTwoSingle.get(0));
            //nextGenerationChild.get(PopulationSize-2).addAll(ElitismTwoSingle.get(1));
            System.out.println("~~~~~~~~~~~迭代次数 "+count+"min fitness~~~ "+min_fitness);
            //w.write(count+":  "+min_fitness+"\n");
            System.out.println("best individual"+nextGenerationChild.get(0));
            /*
            for(int print_index=0;print_index<nextGenerationChild.get(0).size();print_index++)
            {
               // System.out.println((nextGenerationChild.get(0).get(print_index)));
            }
            */
            //System.out.println(nextGenerationChild.get(0));
            //System.out.println(nextGenerationChild.get(PopulationSize/2));
            //System.out.println(nextGenerationChild.get(PopulationSize-1));

        }
        VisitedNumberSet=caculateVistitedNumberSet(nextGenerationChild,VertexALLNumber,VertexPreNumber);
        //对于最终代个体重新计算融合方案的点,边的集合
        VertexNumberAfterFusion=CaculateVertexNumberAfterFusion(nextGenerationChild,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        EdgeNumberAfterFusion=CaculateEdgeNumberAfterFusion(nextGenerationChild,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        fitness=CaculateFitness(nextGenerationChild,VertexALLNumber,Graphset,VertexPreNumber,VisitedNumberSet);
        min_fitness=Long.MAX_VALUE;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            if(fitness[populationIndex]<min_fitness)
            {
                min_fitness=fitness[populationIndex];
                bestIndividualIndex=populationIndex;
            }
        }
        System.out.println("minfitness~~~~~~~~~~~~~~~~~~"+min_fitness);

        //w.write("经过遗传算法优化后的图融合点数"+VertexNumberAfterFusion[bestIndividualIndex]+"\n");
        //w.write("经过遗传算法优化后的图融合边数"+EdgeNumberAfterFusion[bestIndividualIndex]+"\n");
        //w.write("经过遗传算法优化后的最佳个体index: "+bestIndividualIndex+"\n");
        //w.write(generationNumber+" "+min_fitness+"\n");
        //末代最有个体
        w.write(""+min_fitness);
        //System.out.println("best individual final"+nextGenerationChild.get(0));
        //System.out.println("best indivdual final"+nextGenerationChild.get(1));
        //System.out.println("best indivdual final"+nextGenerationChild.get(2));
        //System.out.println("best indivdual final"+nextGenerationChild.get(3));

    }


    //判断新个体与种群中的任意个体是否相似
    public static boolean isExistedInPopulation(List<List<Integer>> single,List<List<List<Integer>>> population)
    {
        for(int i=0;i<PopulationSize;i++)
        {
            if(isSameSingle(single,population.get(i)))
            {
                return true;
            }
        }
        return false;
    }

    //判断种群里面的个体是否相同
    public static boolean isSameSingle(List<List<Integer>>singleA,List<List<Integer>> singleB)
    {
        int row,col;
        for(row=0;row<singleA.size();row++)
        {
            for(col=0;col<singleA.get(row).size();col++)
            {
                if(singleA.get(row).get(col)!=singleB.get(row).get(col))
                {
                    return false;
                }
            }
        }
        return true;
    }
    //再每次迭代的末尾更新,确保相似个体总数不超过DiversityThreshold个
    //keepDiversity操作每次加入新个体时,并没有对该个体的冲突数组进行置换
    public static void keepDiversty(List<List<List<Integer>>> nextGenerationChild,
                                    boolean DiversityVisitedPopulationSet[],int Graphset[][][],int VertexAllNumber,int VertexPreNumber[])
    {
        Random random=new Random();
        //System.out.println("开始保证种群多样性~~~~~~~~~~~");
        int populationIndex=0;
        int populationIndexNext=0;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            List<Integer> count=new ArrayList<>();
            //如果之前已经被访问过
            if(DiversityVisitedPopulationSet[populationIndex])
            {
                continue;
            }
            //对每个个体如果有重复个体
            for(populationIndexNext=0;populationIndexNext<PopulationSize;populationIndexNext++)
            {
                if(isSameSingle(nextGenerationChild.get(populationIndex),nextGenerationChild.get(populationIndexNext)))
                {
                    count.add(populationIndexNext);
                }
            }
            for(Integer key:count)
            {
                DiversityVisitedPopulationSet[key]=true;
            }
            if(count.size()>DiversityThreshold)
            {
                //System.out.println("count.size(): ~~~~~~~~~~~~~"+count.size());
                int excess=count.size()-DiversityThreshold;
                //对超出的重复个体
                for(int i=0;i<excess;i++)
                {
                    //从重复队列里面先一个数组出来
                    int pos=random.nextInt(count.size());
                    List<List<Integer>>  newSingle=generateNewSingle(Graphset,VertexAllNumber,VertexPreNumber);
                    while(isExistedInPopulation(newSingle,nextGenerationChild))
                    {
                        newSingle=generateNewSingle(Graphset,VertexAllNumber,VertexPreNumber);
                    }
                    int replaceSingleNumber=count.get(pos);
                   // System.out.println("由于重复个体太多而替换的个体的索引pos is ~~"+pos);
                    nextGenerationChild.get(replaceSingleNumber).clear();
                    nextGenerationChild.get(replaceSingleNumber).addAll(newSingle);
                }
            }
        }
        //System.out.println("结束保证种群多样性~~~~~~~~~~~");
    }

    public static Comparator<Individual> fitnessComparator=new Comparator<Individual>() {
        @Override
        public int compare(Individual o1, Individual o2) {
            return o1.fitness-o2.fitness;
        }
    };

    //取前两最好的个体,将上一代的parent中最好的两个赋予下一代最后的两个,返回这两个
    //取上一代中最好的top2个体或者top2%的个体,不进行crossOver直接复制到下一代child,但进行变异
    public static Queue<Individual> Elitism( List<List<List<Integer>>> nextGenerationParent,int fitness[],int VertexAllNumber,boolean VisitedNumberSet[][])
    {
       // List<List<List<Integer>>> ElitismSingles=new ArrayList<>();

        Queue<Individual> PopulationPQ=new PriorityQueue<>(PopulationSize,fitnessComparator);
        Queue<Individual> ElitismPQ=new PriorityQueue<>(ElitismNumber,fitnessComparator);

        //构造优先队列,讲所有的个体放入种群优先队列中,重载比较器
        for(int index=0;index<PopulationSize;index++)
        {
            PopulationPQ.offer(new Individual(nextGenerationParent.get(index),fitness[index],VisitedNumberSet[index]));
        }
        //取出topN精英加入精英优先队列
        for(int index=0;index<ElitismNumber;index++)
        {
            //System.out.println("获得当代精英适应度~"+PopulationPQ.peek().fitness);
            ElitismPQ.add(PopulationPQ.poll());
           // ElitismSingles.add(ElitismPQ.poll().graphFusionCode);

        }
        //System.out.println("当代精英种族最精英个体"+ElitismPQ.peek().graphFusionCode);
        return ElitismPQ;


        /*

        int bestfitness=0;
        int second_best_fitness=Integer.MAX_VALUE;
        int bestfitnessIndex=0;
        int second_best_fitness_Index=Integer.MAX_VALUE;
        for(int index=0;index<PopulationSize;index++)
        {
            if(fitness[index]<bestfitness)
            {
                second_best_fitness=bestfitness;
                bestfitness=fitness[index];
                bestfitnessIndex=index;
                continue;
            }
            if(fitness[index]<second_best_fitness)
            {
                second_best_fitness=fitness[index];
                second_best_fitness_Index=index;
                continue;
            }
        }
        for(int i=0;i<PopulationSize/50;i++)
        {
            ElitismTwoSingle.add(nextGenerationParent.get(bestfitnessIndex));
            ElitismTwoSingle.add(nextGenerationParent.get(second_best_fitness_Index));
        }

        return ElitismTwoSingle;
        */
    }


    //精英主义替换,用选择交叉变异后的N个个体加上选择前的E个精英混合,选择其中靠前的N个个体
    public static void ElitismReplace(List<List<List<Integer>>> nextGenerationChild,Queue<Individual> ElitismPQ,
                                      int fitness[],boolean VisitedNumberSet[][],int VertextALLNumber)
    {
        Queue<Individual> ElitismAndPopulationPQ=new PriorityQueue<>(PopulationSize+ElitismNumber,fitnessComparator);
        //将下一代个体全部放入优先队列
        for(int index=0;index<PopulationSize;index++)
        {
            ElitismAndPopulationPQ.offer(new Individual(nextGenerationChild.get(index),fitness[index],VisitedNumberSet[index]));
        }
        //将上一代保存的精英全部放入优先队列
        for(int index=0;index<ElitismNumber;index++)
        {
            //System.out.println("查看上一代保存精英"+ElitismPQ.peek().fitness);
            ElitismAndPopulationPQ.offer(ElitismPQ.poll());
            //ElitismAndPopulationPQ.remove();
        }
        //选择前N个个体进入下一代nextgenerationParents
        for(int index=0;index<PopulationSize;index++)
        {
            nextGenerationChild.get(index).clear();
            fitness[index]=ElitismAndPopulationPQ.peek().fitness;
            //复制冲突数组
            for(int count=0;count<VertextALLNumber;count++)
            {
                VisitedNumberSet[index][count]=ElitismAndPopulationPQ.peek().VisitedNumberSetSingle[count];
            }
            nextGenerationChild.get(index).addAll(ElitismAndPopulationPQ.poll().graphFusionCode);
           // System.out.println("精英主义替换的当代种群编码~~~"+nextGenerationChild.get(index));
            //System.out.println("精英主义替换后的当代种群~~~"+fitness[index]);

        }
        //System.out.println("精英主义替换后的当代种群编码~~~~"+nextGenerationChild);
    }

    //细粒度变异
    public static  void fine_grained_mutation(List<List<List<Integer>>> nextGenerationChild,int Graphset[][][],
                                              int VertexAllNumber,int VertexPreNumber[], boolean VisitedNumberSet[][])
    {
        double mutation_pro_threshhold=Math.random();
        int populationIndex,col;

        int VertexNumberAfterFusion[]=new int[PopulationSize];

        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            //对每个个体,计算其VisitedNumberSet为True的总数
            for(col=0;col<VertexAllNumber;col++)
            {
                //如果为false,说明该列未被置为冲突,VertextNumber++
                if(!VisitedNumberSet[populationIndex][col])
                {
                    VertexNumberAfterFusion[populationIndex]++;
                }
            }

        }
        //System.out.println(nextGenerationChild.get(0));
       // System.out.println("VertextNumberAfterFusion"+VertexNumberAfterFusion[0]);

        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            if(mutation_pro_threshhold<0.15)
            {
                //System.out.println("before mutation"+"发生变异的个体编号是"+populationIndex);
                //System.out.println(nextGenerationChild.get(populationIndex));
                //System.out.println(nextGenerationChild.get(1));
                double mutation_type_choose=Math.random();
                if(mutation_type_choose<1.0/VertexNumberAfterFusion[populationIndex])
                {
                    //System.out.println("发生分裂操作");
                    fine_grained_mutation_split(nextGenerationChild,Graphset,VertexAllNumber,VertexPreNumber,VisitedNumberSet,populationIndex);
                    //System.out.println("分裂操作以后");
                    //System.out.println(nextGenerationChild.get(populationIndex));
                    //System.out.println("After mutation"+"发生变异的个体编号是"+populationIndex+"~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                   // System.out.println(nextGenerationChild.get(1));

                }
                else
                {
                    //System.out.println("发生交换操作");
                    fine_grained_mutation_exchange(nextGenerationChild,Graphset,VertexAllNumber,VertexPreNumber,VisitedNumberSet,populationIndex);
                    //System.out.println("交换操作以后");
                    //System.out.println(nextGenerationChild.get(populationIndex));
                    //System.out.println("After mutation"+"发生变异的个体编号是"+populationIndex+"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");

                }
            }


        }



        return ;
    }


    //判断融合编码的某一列长度是否大于1
    public static int caculateFusionCodeColLength(List<List<Integer>> nextGenerationChildSingle,int choose_mutation_col,int VertexPreNumber[])
    {
        int belong_graph=belongtoWhichGraph(choose_mutation_col,VertexPreNumber);
        //该列一共有多长
        int choose_mutation_col_height_all=GraphNumber-1-belong_graph;
        //定义有意义的(非-1的)该列长度,默认为1
        int valid_code_col_length=1;
        //从belong_graph开始遍历该列
        for(int row=belong_graph;row<GraphNumber-1;row++)
        {
            if(nextGenerationChildSingle.get(row).get(choose_mutation_col)!=-1)
            {
                valid_code_col_length++;
            }
        }

        return valid_code_col_length;
    }

    //细粒度变异之融合策略
    //从
    public static void fine_grained_mutation_fusion(List<List<List<Integer>>> nextGenerationChild,int Graphset[][][],
                                                    int VertexAllNumber,int VertexPreNumber[],boolean VisitedNumberSet[][], int populationIndex)
    {
        int choose_mutation_col;
        List<List<Integer>> tempGenerationChildSingle=nextGenerationChild.get(populationIndex);
        Random random=new Random();
        //选择要融合的某一列
        choose_mutation_col=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
    }

    //细粒度变异之融合
    /*
    public static void fine_grained_mutation_fusion(List<List<List<Integer>>> nextGenerationChild,int Graphset[][][],
                                                    int VertexAllNumber,int VertexPreNumber[],boolean VisitedNumberSet[][]) {
        int populationIndex;
        //每次循环内部参数申明
        int choose_mutation_col, belong_graph, choose_mutation_col_height_all, choose_mutation_col_height, choose_mutation_row;
        for (populationIndex = 0; populationIndex < PopulationSize; populationIndex++) {
            List<List<Integer>> tempGenerationChildSingle = nextGenerationChild.get(populationIndex);
            Random random = new Random();
            //选中某一列变异
            choose_mutation_col = random.nextInt(tempGenerationChildSingle.get(GraphNumber - 2).size());
            while (VisitedNumberSet[populationIndex][choose_mutation_col]) {
                choose_mutation_col = random.nextInt(tempGenerationChildSingle.get(GraphNumber - 2).size());
            }
            belong_graph = belongtoWhichGraph(choose_mutation_col, VertexPreNumber);
            choose_mutation_col_height_all = GraphNumber - 1 - belong_graph;
            //变异列从起始行到行高
            choose_mutation_col_height = random.nextInt(choose_mutation_col_height_all);
            //实际最终行数
            choose_mutation_row = belong_graph + choose_mutation_col_height;
            //  System.out.println("belong_graph"+belong_graph);
            // System.out.println("choose_mutation_col "+choose_mutation_col);
            //System.out.println("choose_mutation_row_height"+choose_mutation_col_height);
            //System.out.println("cchoose_mutation_row" +choose_mutation_row);

            //如果选中的该点不为-1,则分裂
            if (tempGenerationChildSingle.get(choose_mutation_row).get(choose_mutation_col) != -1) {
                int swap_col = VertexPreNumber[choose_mutation_row + 1] + tempGenerationChildSingle.get(choose_mutation_row).get(choose_mutation_col);
                VisitedNumberSet[populationIndex][swap_col] = false;
                tempGenerationChildSingle.get(choose_mutation_row).set(choose_mutation_col, -1);
                //  System.out.println("mutation VisitedNumberSet~~~");
            }
        }
    }
    */


    //细粒度变异之分裂
    public static void fine_grained_mutation_split(List<List<List<Integer>>> nextGenerationChild,int Graphset[][][],
                                                   int VertexAllNumber,int VertexPreNumber[],boolean VisitedNumberSet[][], int populationIndex)
    {

        //每次循环内部参数申明
        int choose_mutation_col,belong_graph,choose_mutation_col_height_all,choose_mutation_col_height,choose_mutation_row;
            List<List<Integer>> tempGenerationChildSingle=nextGenerationChild.get(populationIndex);
            Random random=new Random();
            choose_mutation_col=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
            while(VisitedNumberSet[populationIndex][choose_mutation_col])
            {
                choose_mutation_col=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
            }
            belong_graph=belongtoWhichGraph(choose_mutation_col,VertexPreNumber);
            //choose_mutation_col_height_all=GraphNumber-1-belong_graph;
            //变异列从起始行到行高
            //choose_mutation_col_height=random.nextInt(choose_mutation_col_height_all);

            //System.out.println("belong_graph"+belong_graph);
            //System.out.println("choose_mutation_col "+choose_mutation_col);
            //System.out.println("choose_mutation_row_height"+choose_mutation_col_height);
            //System.out.println("cchoose_mutation_row" +choose_mutation_row);
            int choose_mutation_col_length=caculateFusionCodeColLength(tempGenerationChildSingle,choose_mutation_col,VertexPreNumber);

            //如果选中的列长度(即融合后的集合长度<=1),则该列无法分裂,返回
            if(choose_mutation_col_length<=1)
            {
                //System.out.println("采用分裂方式1,基数为1的集合无法分裂");
                return;
            }
           //选出所有选中列中有有效元素的row
            List<Integer> valid_code_list=new ArrayList<>();



            for(int valid_row=belong_graph;valid_row<GraphNumber-1;valid_row++)
            {
                if(tempGenerationChildSingle.get(valid_row).get(choose_mutation_col)!=-1)
                {
                    valid_code_list.add(valid_row);
                }
            }

            //该集合中最后有-1(即该列对应的图列首元素)
            valid_code_list.add(-1);

            //System.out.println(valid_code_list);


           //定义选择分裂的那一行的number
            int choose_mutation_row_index=random.nextInt(valid_code_list.size());
            int choose_mutation_row_number=valid_code_list.get(choose_mutation_row_index);

           //实际最终行数
           //choose_mutation_row=belong_graph+choose_mutation_col_height;
          //如果number=-1,则将该中其他元素全部分离出去
            if(choose_mutation_row_number==-1)
            {
                //System.out.println("采用分裂方式二,保留该列首位元素,将其他列元素迁移到第一个元素所在的列");
                //System.out.println(tempGenerationChildSingle);
                //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                int choose_mutation_row_first=valid_code_list.get(0);
                int swap_col=VertexPreNumber[choose_mutation_row_first+1]+tempGenerationChildSingle.get(choose_mutation_row_first).get(choose_mutation_col);
                tempGenerationChildSingle.get(choose_mutation_row_first).set(choose_mutation_col,-1);
                //将list里面的元素全部取出来与对应的swap_col进行交换
                 for(int index=1;index<=valid_code_list.size()-2;index++)
                 {
                     int choose_mutation_row_second=valid_code_list.get(index);
                     int choose_mutation_split_number=tempGenerationChildSingle.get(choose_mutation_row_second).get(choose_mutation_col);
                     tempGenerationChildSingle.get(choose_mutation_row_second).set(swap_col,choose_mutation_split_number);
                     //将原来列的值置为-1
                     tempGenerationChildSingle.get(choose_mutation_row_second).set(choose_mutation_col,-1);
                 }
                VisitedNumberSet[populationIndex][swap_col]=false;
                //VisitedNumberSet[populationIndex][choose_mutation_col]=true;
               // System.out.println("分裂方式二以后");
                //System.out.println(tempGenerationChildSingle);
                //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));

            }
            //如果number不为-1,则将该列对应的元素分裂出去
            else
            {
                //System.out.println("采用分裂方式三,将该列不为首位的元素分裂出去,即将该元素对应的列置为Visited=false");
                //System.out.println(tempGenerationChildSingle);
                //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
               // System.out.println(choose_mutation_row_number+"row_number");
                //System.out.println(choose_mutation_col+"col");
                int swap_col=VertexPreNumber[choose_mutation_row_number+1]+tempGenerationChildSingle.get(choose_mutation_row_number).get(choose_mutation_col);
                VisitedNumberSet[populationIndex][swap_col]=false;
                tempGenerationChildSingle.get(choose_mutation_row_number).set(choose_mutation_col,-1);
                //System.out.println("分裂方式三以后~~~~~~~~~~~~~~");
                //System.out.println(tempGenerationChildSingle);
                //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
            }

            //如果选中的该点不为-1,则分裂
        /*
            if(tempGenerationChildSingle.get(choose_mutation_row).get(choose_mutation_col)!=-1)
            {
                int swap_col=VertexPreNumber[choose_mutation_row+1]+tempGenerationChildSingle.get(choose_mutation_row).get(choose_mutation_col);
                VisitedNumberSet[populationIndex][swap_col]=false;
                tempGenerationChildSingle.get(choose_mutation_row).set(choose_mutation_col,-1);
                //System.out.println("mutation VisitedNumberSet~~~");
            }
         */
        return;
    }
    //细粒度变异之交换
    //交换必须建立在两个集合的某一元素都是有效元素的前提下
    public static void fine_grained_mutation_exchange(List<List<List<Integer>>> nextGenerationChild,int Graphset[][][],
                                                      int VertexAllNumber,int VertexPreNumber[],boolean VisitedNumberSet[][],int populationIndex)
    {

            List<List<Integer>> tempGenerationChildSingle=nextGenerationChild.get(populationIndex);
            Random random=new Random();
            int col_swap_left=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
            //确保left列不为冲突列
            while(VisitedNumberSet[populationIndex][col_swap_left])
            {
                col_swap_left=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
            }

            //确保right列不为冲突列
            int col_swap_right=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
            while(VisitedNumberSet[populationIndex][col_swap_right]||(col_swap_left==col_swap_right))
            {
                col_swap_right=random.nextInt(tempGenerationChildSingle.get(GraphNumber-2).size());
            }

            //如果left>right,left与right交换
            if(col_swap_left>col_swap_right)
            {
                int temp=col_swap_left;
                col_swap_left=col_swap_right;
                col_swap_right=temp;
            }


        //判断右边这一列对应的元素有效编码的长度
        int col_swap_right_code_length=caculateFusionCodeColLength(tempGenerationChildSingle,col_swap_right,VertexPreNumber);

         /*
            //如果右边位于第一个图
           if(col_swap_right<Graphset[0].length)
           {

                   //int map_graph_index_left=col_swap_right;
                   //int map_graph_index_right=col_swap_left;
                   for(int row=0;row<GraphNumber-1;row++)
                   {
                       int swaptemp=tempGenerationChildSingle.get(row).get(col_swap_left);
                       tempGenerationChildSingle.get(row).set(col_swap_left,tempGenerationChildSingle.get(row).get(col_swap_right));
                       tempGenerationChildSingle.get(row).set(col_swap_right,swaptemp);
                   }

               return ;
           }
           */
           int belong_graph_left=belongtoWhichGraph(col_swap_left,VertexPreNumber);
           int belong_graph_right=belongtoWhichGraph(col_swap_right,VertexPreNumber);
           //如果右边该列对应的集合只有一个节点
            if(col_swap_right_code_length==1)
            {

                //对应的left的映射行
                int map_left_row=belong_graph_right-1;
                //如果右列和左列属于同一个图(编码同一行),直接对下方编码交换
                if(belong_graph_left==belong_graph_right)
                {
                   // System.out.println("变异交换策略方式1,左右列队首位于同一个图,右列只有队首元素,则将左列元素全部交换到右列即可");
                    //System.out.println(tempGenerationChildSingle);
                    //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                    for(int belong_graph_row_index=belong_graph_left;belong_graph_row_index<GraphNumber-1;belong_graph_row_index++)
                    {
                        int swaptemp=tempGenerationChildSingle.get(belong_graph_row_index).get(col_swap_left);
                        tempGenerationChildSingle.get(belong_graph_row_index).set(col_swap_left,tempGenerationChildSingle.get(belong_graph_row_index).get(col_swap_right));
                        tempGenerationChildSingle.get(belong_graph_row_index).set(col_swap_right,swaptemp);
                    }
                    //System.out.println(tempGenerationChildSingle);
                    //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                    //System.out.println("变异方式1结束");
                }
                //如果右列和左列不属于同一个图(即左列位于右列的图前方)
                else
                {


                    //如果对应left的映射编码为-1(左集合内没有该图对应的节点),将该列置为Visited True,将left列编码设为对应
                   // System.out.println("map_left_row: "+map_left_row);
                    //System.out.println("belong_graph"+belong_graph_right);
                    //System.out.println("col_swap_left"+col_swap_left);
                    if(tempGenerationChildSingle.get(map_left_row).get(col_swap_left)==-1)
                    {
                        //System.out.println("变异交换策略方式2,左右列队首没有位于同一个图,右列只有队首元素,左列无对应元素" +
                        //        "则将右列元素置为true,右列元素来到左列,左列对应元素置为false");
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        int map_graph_index=col_swap_right-VertexPreNumber[belong_graph_right];
                        tempGenerationChildSingle.get(map_left_row).set(col_swap_left,map_graph_index);
                        VisitedNumberSet[populationIndex][col_swap_right]=true;
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        //System.out.println("变异方式2结束");

                    }
                    //如果对应left的映射编码不为-1(左集合内有该图对应的节点),
                    else
                    {
                        //System.out.println("变异交换策略方式3,左右列队首没有位于同一个图,右列只有队首元素,左列有对应元素" +
                        //"则将右列元素置为true,右列元素来到左列,左列对应元素置为false");
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        ////将左集合的该列置为Visited false,右集合的原列置为Visited True
                        int map_graph_index_left=VertexPreNumber[belong_graph_right]+tempGenerationChildSingle.get(map_left_row).get(col_swap_left);
                        int map_graph_index_right=col_swap_right-VertexPreNumber[belong_graph_right];
                        tempGenerationChildSingle.get(map_left_row).set(col_swap_left,map_graph_index_right);
                        VisitedNumberSet[populationIndex][col_swap_right]=true;
                        VisitedNumberSet[populationIndex][map_graph_index_left]=false;
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        //System.out.println("变异方式3结束");

                    }
                }


                return;
            }



            //如果右边集合的大小|S|>1
            //对右边的那列选出某一个有效的编码作为交换
            //int belong_graph=belongtoWhichGraph(col_swap_right,VertexPreNumber);
            int map_right_head_row=belong_graph_right-1;
            //选出所有选中右列中有有效元素的row
            List<Integer> valid_code_list=new ArrayList<>();
            for(int valid_row=belong_graph_right;valid_row<GraphNumber-1;valid_row++)
            {
                if(tempGenerationChildSingle.get(valid_row).get(col_swap_right)!=-1)
                {
                    valid_code_list.add(valid_row);
                }
            }
            //该集合中最后有-1(即该列对应的图列首元素)
            valid_code_list.add(-1);

            //定义选择分裂的那一行的number
            int choose_mutation_row_index=random.nextInt(valid_code_list.size());
            int choose_mutation_row_number=valid_code_list.get(choose_mutation_row_index);
            //如果选中了头部元素
            if(choose_mutation_row_number==-1)
            {

                //如果左右列位于同一个图,则交换下方编码
                if(belong_graph_left==belong_graph_right)
                {
                   // System.out.println("变异交换策略方式4,左右列队首位于同一个图,右列长度大于1,选中右列头部" +
                            //"则将左右列元素互换即可");
                   // System.out.println(tempGenerationChildSingle);
                    //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                    for(int belong_graph_row_index=belong_graph_left;belong_graph_row_index<GraphNumber-1;belong_graph_row_index++)
                    {
                        int swaptemp=tempGenerationChildSingle.get(belong_graph_row_index).get(col_swap_left);
                        tempGenerationChildSingle.get(belong_graph_row_index).set(col_swap_left,tempGenerationChildSingle.get(belong_graph_row_index).get(col_swap_right));
                        tempGenerationChildSingle.get(belong_graph_row_index).set(col_swap_right,swaptemp);
                    }
                    //System.out.println(tempGenerationChildSingle);
                    //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                    //System.out.println("变异方式4结束");
                }

                //如果左右列不位于同一个图
                else
                {

                    //如果对应left的映射编码为-1
                    if(tempGenerationChildSingle.get(map_right_head_row).get(col_swap_left)==-1)
                    {
                        //System.out.println("变异交换策略方式5,左右列队首没有位于同一个图,右列长度大于1,选中右列头部" +
                                //"左列没有该图对应元素,则将右列元素移过来,右列其他元素移动到右列第一个元素对应的列");
                        //将除头部以外的有效code交换到右半部分
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        //取出第一个元素,找到它对应的VisitedNumberset
                        int choose_mutation_row_first=valid_code_list.get(0);
                        int swap_col_after_right=VertexPreNumber[choose_mutation_row_first+1]+tempGenerationChildSingle.get(choose_mutation_row_first).get(col_swap_right);
                        //将第一个元素置为-1
                        tempGenerationChildSingle.get(choose_mutation_row_first).set(col_swap_right,-1);
                        //将list里面的元素全部取出来与对应的swap_col进行交换
                        for(int index=1;index<=valid_code_list.size()-2;index++)
                        {
                            int choose_mutation_row_second=valid_code_list.get(index);
                            int choose_mutation_exchange_number=tempGenerationChildSingle.get(choose_mutation_row_second).get(col_swap_right);
                            tempGenerationChildSingle.get(choose_mutation_row_second).set(swap_col_after_right,choose_mutation_exchange_number);
                            //将原来列的值置为-1
                            tempGenerationChildSingle.get(choose_mutation_row_second).set(col_swap_right,-1);
                        }
                        //将right列置为已访问,将right列第一个不为-1的数对应的列置为false未访问
                        VisitedNumberSet[populationIndex][col_swap_right]=true;
                        VisitedNumberSet[populationIndex][swap_col_after_right]=false;
                        //进行交换
                        tempGenerationChildSingle.get(map_right_head_row).set(col_swap_left,col_swap_right-VertexPreNumber[belong_graph_right]);
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        //System.out.println("变异方式5结束");
                    }
                    //如果对应left的映射编码不为-1
                    else
                    {
                        //System.out.println("变异交换策略方式6,左右列队首没有位于同一个图,右列长度大于1,选中右列头部" +
                                //"左列有该图对应元素,则将右列元素移过来,右列其他元素移动到左列该元素对应的列");
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        int map_left_col=VertexPreNumber[belong_graph_right]+tempGenerationChildSingle.get(belong_graph_right-1).get(col_swap_left);
                        //如果两边都有有效编码,将右列的有效编码交换到对应的左列映射编码,然后将右列头部填入左列,修改冲突数组
                        for(int row_index=belong_graph_right;row_index<GraphNumber-1;row_index++)
                        {
                            //int swaptemp=tempGenerationChildSingle.get(row_index).get(col_swap_left);
                            tempGenerationChildSingle.get(row_index).set(map_left_col,tempGenerationChildSingle.get(row_index).get(col_swap_right));
                            tempGenerationChildSingle.get(row_index).set(col_swap_right,-1);
                            //tempGenerationChildSingle.get(row_index).set(col_swap_right,swaptemp);
                        }
                        //最后将右列的队首元素值填入左列对应行
                        tempGenerationChildSingle.get(map_right_head_row).set(col_swap_left,col_swap_right-VertexPreNumber[belong_graph_right]);
                        VisitedNumberSet[populationIndex][col_swap_right]=true;
                        VisitedNumberSet[populationIndex][map_left_col]=false;
                        //System.out.println(tempGenerationChildSingle);
                        //System.out.println(Arrays.toString(VisitedNumberSet[populationIndex]));
                        //System.out.println("变异方式6结束");

                    }

                }
            }

            //如果没有选中头部元素
            else
            {
               // System.out.println("变异交换策略方式7,左右列队首没有位于同一个图,右列长度大于1,未选中右列头部" +
                 //       "则直接对左右列元素进行交换即可");
                //System.out.println(tempGenerationChildSingle);
                int swaptemp=tempGenerationChildSingle.get(choose_mutation_row_number).get(col_swap_left);
                tempGenerationChildSingle.get(choose_mutation_row_number).set(col_swap_left,tempGenerationChildSingle.get(choose_mutation_row_number).get(col_swap_right));
                tempGenerationChildSingle.get(choose_mutation_row_number).set(col_swap_right,swaptemp);
                //System.out.println(tempGenerationChildSingle);

            }

            /*
            int col_right_heigh_all=GraphNumber-1-belongtoWhichGraph(col_swap_right,VertexPreNumber);
            int col_right_height=random.nextInt(col_right_heigh_all);
            int swap_row=belongtoWhichGraph(col_swap_right,VertexPreNumber)+col_right_height;
            */
            //细粒度交换,在swap_row 行对col_right和col_left列做交换(并没有判断是否交换的编码为-1)


        return;

    }

    //粗粒度交叉
    //先交换编码和对应的VisitedNumberSet,
    public static  List<List<List<Integer>>> coarse_grained_crossover(List<List<List<Integer>>> nextGenerationParent,
                                                                      int VertexPreNumber[],boolean VisitedNumberSet[][],int VertexALLNumber)
    {
        List<List<List<Integer>>> nextGenerationChild=new ArrayList<>();
        int populationIndex;
        for (populationIndex=0;populationIndex<PopulationSize;populationIndex+=2)
        {
            double cross_over_chosse_prob=Math.random();
            //生成子代1,2
            List<List<Integer>> nextGenerationChildSingleLeft=new ArrayList<>();
            List<List<Integer>> nextGenerationChildSingleRight=new ArrayList<>();



            //左上+右下结合,右上和左下结合,交叉,发生概率较高
            if(cross_over_chosse_prob<1.0)
            {
                for(int row_up=0;row_up<(GraphNumber-1)/2;row_up++)
                {
                    nextGenerationChildSingleLeft.add(nextGenerationParent.get(populationIndex).get(row_up));
                    nextGenerationChildSingleRight.add(nextGenerationParent.get(populationIndex+1).get(row_up));
                }
                for(int row_down=(GraphNumber-1)/2;row_down<GraphNumber-1;row_down++)
                {
                    nextGenerationChildSingleLeft.add(nextGenerationParent.get(populationIndex+1).get(row_down));
                    nextGenerationChildSingleRight.add(nextGenerationParent.get(populationIndex).get(row_down));
                }
                int row_up_end=(GraphNumber-1)/2;

                //对应的图+1,对Visted已访问列数组进行交换
                int swap_split_col=VertexPreNumber[row_up_end+1];

                //System.out.println("modify before VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
                //System.out.println("modify before VisitedSet"+Arrays.toString(VisitedNumberSet[1]));
                //System.out.println("swap split col is~~~~~~~~~~"+swap_split_col);
                for(int col_index=swap_split_col;col_index<VertexALLNumber;col_index++)
                {
                    boolean tempBool=VisitedNumberSet[populationIndex][col_index];
                    VisitedNumberSet[populationIndex][col_index]=VisitedNumberSet[populationIndex+1][col_index];
                    VisitedNumberSet[populationIndex+1][col_index]=tempBool;
                }
                //System.out.println("modify after VisitedSet"+Arrays.toString(VisitedNumberSet[0]));
                //System.out.println("modify after VisitedSet"+Arrays.toString(VisitedNumberSet[1]));
                //System.out.println("modify before ~~~~~~`");
                //System.out.println(nextGenerationChildSingleLeft);
                //System.out.println(nextGenerationChildSingleRight);
                ModifyNextGenerationChildSingle(nextGenerationChildSingleLeft,VisitedNumberSet,row_up_end,swap_split_col,VertexPreNumber,populationIndex);
                ModifyNextGenerationChildSingle(nextGenerationChildSingleRight,VisitedNumberSet,row_up_end,swap_split_col,VertexPreNumber,populationIndex);
                //System.out.println("modify After ~~~~~");
                //System.out.println(nextGenerationChildSingleLeft);
                //System.out.println(nextGenerationChildSingleRight);
            }
            //左上加左上结合,右上和右下结合,类似于不交叉
            else
            {
                nextGenerationChildSingleLeft.addAll(nextGenerationParent.get(populationIndex));
                nextGenerationChildSingleRight.addAll(nextGenerationParent.get(populationIndex+1));
            }
            int middle_up=(GraphNumber-1)/2-1;
            //System.out.print("modify before ");
            //System.out.println(nextGenerationChildSingleLeft);
            // ModifyNextGenerationChildSingle(nextGenerationChildSingleLeft,middle_up,middle_up+1,VertexPreNumber);
            //ModifyNextGenerationChildSingle(nextGenerationChildSingleLeft,middle_up,middle_up+1,VertexPreNumber);
            //System.out.print("modify After ");
            //System.out.println(nextGenerationChildSingleLeft);
            //System.out.println(nextGenerationChildSingleRight);
            nextGenerationChild.add(nextGenerationChildSingleLeft);
            nextGenerationChild.add(nextGenerationChildSingleRight);
        }
        //System.out.println("Afeter CrossOver ~~~~~~~~~~~~~~~~~~~~");
        //System.out.println(nextGenerationChild.get(0));
        //System.out.println(nextGenerationChild.get(1));
        return nextGenerationChild;
    }


    //粗粒度交叉后修正
    public static void ModifyNextGenerationChildSingle(List<List<Integer>> nextGenerationChildSingle,boolean VisitedNumberSet[][],
                                                       int middle_up_row,int split_col_right_begin,int VertexPreNumber[],
                                                       int populationIndex)
    {
        int  row,col;
        int belong_row=middle_up_row;
        //对块交换后的第一行进行判定,如果第一行没有违背,则下面的行都不会违背
        for(col=0;col<nextGenerationChildSingle.get(belong_row).size();col++)
        {
            //如果列Visted为True,需找到第列第一个不为-1的,进行交换
            if(VisitedNumberSet[populationIndex][col])
            {
                int swap_target_col=-1;
                for(row=belong_row;row<GraphNumber-1;row++)
                {
                    //如果没到最后一行&&没找到对应列
                    if(row!=GraphNumber-2&&swap_target_col<0)
                    {
                        if((nextGenerationChildSingle.get(row).get(col)==-1))
                        {
                            // System.out.println("swap11`~~~~");
                            continue;
                        }
                        //如果在该行找到不为-1的元素,则记录下该元素所对应列
                        else if(nextGenerationChildSingle.get(row).get(col)!=-1)
                        {
                            //System.out.println("swap12`~~~~");
                            swap_target_col=VertexPreNumber[row+1]+nextGenerationChildSingle.get(row).get(col);
                        }

                    }
                    //如果到了最后一行
                    else if(row==GraphNumber-2)
                    {
                        if(nextGenerationChildSingle.get(row).get(col)!=-1)
                        {
                            swap_target_col=VertexPreNumber[row+1]+nextGenerationChildSingle.get(row).get(col);
                            //修正已访问数组
                            nextGenerationChildSingle.get(row).set(col,-1);
                            VisitedNumberSet[populationIndex][swap_target_col]=false;
                            // System.out.println("swap21`~~~~");
                        }
                        // System.out.println("swap2`~~~~");
                    }
                    //如果中途有不为-1的,则与对应的列完全交换,原列全置为-1,内部执行完以后直接break;
                    else if(swap_target_col>=0)
                    {
                        for(int row_index_2=row-1;row_index_2<GraphNumber-2;row_index_2++)
                        {
                            nextGenerationChildSingle.get(row_index_2+1).set(swap_target_col,nextGenerationChildSingle.get(row_index_2).get(col));
                            nextGenerationChildSingle.get(row_index_2).set(col,-1);
                        }
                        //最后一行同样置为-1
                        nextGenerationChildSingle.get(GraphNumber-2).set(col,-1);
                        VisitedNumberSet[populationIndex][swap_target_col]=false;
                        //System.out.println("swap3~~~~~");
                        break;
                    }

                }
            }
        }

    }

    //计算VisitedNumberSet数组
    public static boolean[][] caculateVistitedNumberSet(List<List<List<Integer>>> nextGenerationParent,int VertextALLNumber,
                                                        int VertexPreNumber[])
    {
        boolean VisitedNumberSet[][]=new boolean[PopulationSize][VertextALLNumber];
        int i,j,col,row;
        //数组初始化
        for(i=0;i<PopulationSize;i++)
        {
            for(j=0;j<VertextALLNumber;j++)
            {
                VisitedNumberSet[i][j]=false;
            }
        }
        int populationIndex;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            //对单个图融合方案编码
            List<List<Integer>> nextGenerationParentSingle=nextGenerationParent.get(populationIndex);
            for(row=0;row<GraphNumber-1;row++)
            {
                for(col=0;col<nextGenerationParentSingle.get(row).size();col++)
                {
                    //如果该值不为-1,将其对应的列置为Visited
                    if(nextGenerationParentSingle.get(row).get(col)!=-1)
                    {
                        int belong_col=VertexPreNumber[row+1]+nextGenerationParentSingle.get(row).get(col);
                        VisitedNumberSet[populationIndex][belong_col]=true;
                    }
                }
            }
        }
        return VisitedNumberSet;
    }


    ////采用round-wheel轮盘发选出下一代parent群体
    public static List<List<List<Integer>>> Select(int fitness[],List<List<List<Integer>>> firstPopulation)
    {
        List<List<List<Integer>>> nextGenerationParent=new ArrayList<>();
        double prob_fitness_reverse[]=new double[PopulationSize];
        double prob_fitness_reverse_sum=0.0;
        Random random=new Random();
        //盘子比例
        double wheel[]=new double[PopulationSize];
        int populationIndex;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            prob_fitness_reverse[populationIndex]=1.0/fitness[populationIndex];
        }
        //计算fitness求倒后的总概率和
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            prob_fitness_reverse_sum+=prob_fitness_reverse[populationIndex];
        }
        //修正prob_fitness,使之归一化
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            prob_fitness_reverse[populationIndex]=prob_fitness_reverse[populationIndex]/prob_fitness_reverse_sum;
        }

        //计算每个个体占比(选中概率)
        wheel[0]=prob_fitness_reverse[0];
        for(populationIndex=0;populationIndex<PopulationSize-1;populationIndex++)
        {
            wheel[populationIndex+1]=wheel[populationIndex]+prob_fitness_reverse[populationIndex+1];
        }
        // System.out.println("wheels: "+Arrays.toString(wheel));
        //按照概率随机选择填充生成PopulationSize个下一代parents
        //System.out.println("last wheel~~~"+wheel[PopulationSize-1]);
        for(int count=0;count<PopulationSize;count++)
        {
            double ParentProb=Math.random();
            List<List<Integer>> nextGenerationParentSingle=new ArrayList<>();
            //如果生成数稍微越界,则选择最后的那个个体
            if(ParentProb>wheel[PopulationSize-1])
            {
                nextGenerationParentSingle.addAll(firstPopulation.get(PopulationSize-1));
                nextGenerationParent.add(nextGenerationParentSingle);
            }

            for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
            {
                //如果落入该轮盘
                if(ParentProb<wheel[populationIndex])
                {
                    nextGenerationParentSingle.addAll(firstPopulation.get(populationIndex));
                    break;
                }
            }
            nextGenerationParent.add(nextGenerationParentSingle);
        }


        return  nextGenerationParent;
    }

    //判断某列属于哪一行(哪个图),从0开始
    public static int belongtoWhichGraph(int col,int VertexPreNumber[])
    {
        int k=0;
        while(k<GraphNumber&&col>=VertexPreNumber[k])
        {
            k++;
        }
        return k-1;
    }

    //计算融合后的种群适应度
    public static int[] CaculateVertexNumberAfterFusion(List<List<List<Integer>>> firstPopulation,int VertexAllNumber,int Graphset[][][],
                                                        int VertexPreNumber[],boolean VisitedNumbersSet[][])
    {
        int VertexNumberAfterFusion[]=new int[PopulationSize];
        Arrays.fill(VertexNumberAfterFusion,0);
        int populationIndex,row,col,graphIndex,i,j,x,y;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            //对每个个体,计算其VisitedNumberSet为True的总数
            for(col=0;col<VertexAllNumber;col++)
            {
                //如果为false,说明该列未被置为冲突,VertextNumber++
                if(!VisitedNumbersSet[populationIndex][col])
                {
                    VertexNumberAfterFusion[populationIndex]++;
                }
            }

        }
        return VertexNumberAfterFusion;
    }



    //用广义熵作为适应度函数
    public static int[] caculateFitnessEntropy(List<List<List<Integer>>> firstPopulation,int VertexAllNumber,int Graphset[][][],
                                                      int VertexPreNumber[],boolean VisitedNumbersSet[][]) {
        int fitness[] = new int[PopulationSize];
        //融合后的方案图,用邻接矩阵的方式存储,每个矩阵
        Set<Integer> GraphAfterFusion[][][] = new Set[PopulationSize][VertexAllNumber][VertexAllNumber];
        int populationIndex, row, col, graphIndex, i, j, x, y;
        for (populationIndex = 0; populationIndex < PopulationSize; populationIndex++) {
            for (graphIndex = 0; graphIndex < GraphNumber; graphIndex++) {
                for (i = 0; i < Graphset[graphIndex].length; i++) {
                    for (j = i; j < Graphset[graphIndex].length; j++) {
                        //如果图集合原图有边,对于出边和入边点i,j分别找到在融合方案中对应的index
                        if (Graphset[graphIndex][i][j] == 1) {
                            int i_index = 0;
                            int j_index = 0;
                            if (graphIndex == 0) {
                                GraphAfterFusion[populationIndex][i][j].add(0);
                                GraphAfterFusion[populationIndex][j][i].add(0);
                                i_index = i;
                                j_index = j;

                            } else if (graphIndex >= 1) {

                                //对于(graphindex-1)list,如果存在对应关系
                                if (firstPopulation.get(populationIndex).get(graphIndex - 1).contains(i)) {
                                    i_index = firstPopulation.get(populationIndex).get(graphIndex - 1).indexOf(i);

                                }
                                //对于graph-1 list,没有对应关系,则该节点为一列的开头节点
                                else {
                                    i_index = VertexPreNumber[graphIndex] + i;
                                }

                                //对于j同理
                                if (firstPopulation.get(populationIndex).get(graphIndex - 1).contains(j)) {
                                    j_index = firstPopulation.get(populationIndex).get(graphIndex - 1).indexOf(j);

                                }
                                //对于graph-1 list,没有对应关系,则该节点为一列的开头节点
                                else {
                                    j_index = VertexPreNumber[graphIndex] + j;
                                }
                                GraphAfterFusion[populationIndex][i_index][j_index].add(graphIndex);
                                GraphAfterFusion[populationIndex][j_index][i_index].add(graphIndex);

                            }
                            //System.out.println("i_index is: "+i_index +" "+"j_index is"+j_index+" populationIndexis :"+populationIndex);
                        }
                    }
                }
            }

        }
        return fitness;
        //得到元素为集合的邻接矩阵,对每一个融合节点计算熵,先计算每个节点的出边的集合构成,即除当前点外其他点与该点的边连接情况

    }

    public static int[] CaculateEdgeNumberAfterFusion(List<List<List<Integer>>> firstPopulation,int VertexAllNumber,int Graphset[][][],
                                                      int VertexPreNumber[],boolean VisitedNumbersSet[][])
    {
        int EdgeNumberAfterFusion[]=new int[PopulationSize];
        int populationIndex,row,col,graphIndex,i,j,x,y;
        //计算该种群每个编码对应的边的个数,采用邻接矩阵编码
        int GraphAfterFusion[][][]=new int[PopulationSize][VertexAllNumber][VertexAllNumber];
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            for(graphIndex=0;graphIndex<GraphNumber;graphIndex++)
            {
                for(i=0;i<Graphset[graphIndex].length;i++)
                {
                    for(j=i;j<Graphset[graphIndex].length;j++)
                    {
                        //如果图集合原图有边,对于i,j分别修正融合后的方案图矩阵
                        if(Graphset[graphIndex][i][j]==1)
                        {
                            int i_index=0;
                            int j_index=0;
                            if(graphIndex==0)
                            {
                                GraphAfterFusion[populationIndex][i][j]=1;
                                GraphAfterFusion[populationIndex][j][i]=1;
                                i_index=i;
                                j_index=j;

                            }
                            else if(graphIndex>=1)
                            {

                                //对于(graphindex-1)list,如果存在对应关系
                                if(firstPopulation.get(populationIndex).get(graphIndex-1).contains(i))
                                {
                                    i_index=firstPopulation.get(populationIndex).get(graphIndex-1).indexOf(i);

                                }
                                //对于graph-1 list,没有对应关系,则该节点为一列的开头节点
                                else
                                {
                                    i_index=VertexPreNumber[graphIndex]+i;
                                }

                                //对于j同理
                                if(firstPopulation.get(populationIndex).get(graphIndex-1).contains(j))
                                {
                                    j_index=firstPopulation.get(populationIndex).get(graphIndex-1).indexOf(j);

                                }
                                //对于graph-1 list,没有对应关系,则该节点为一列的开头节点
                                else
                                {
                                    j_index=VertexPreNumber[graphIndex]+j;
                                }
                                GraphAfterFusion[populationIndex][i_index][j_index]=1;
                                GraphAfterFusion[populationIndex][j_index][i_index]=1;

                            }
                            //System.out.println("i_index is: "+i_index +" "+"j_index is"+j_index+" populationIndexis :"+populationIndex);
                        }
                    }
                }
            }
            for(x=0;x<VertexAllNumber;x++)
            {
                for(y=x;y<VertexAllNumber;y++)
                {
                    if(GraphAfterFusion[populationIndex][x][y]==1)
                    {
                        EdgeNumberAfterFusion[populationIndex]++;
                    }
                }
            }
        }
        return EdgeNumberAfterFusion;
    }




    //对于种群的适应度,采用|V|*|E|作为计算方式,需要种群编码,原图集合作为输入
    //fitness的计算取决与VistedNumberSet是否正确
    public static int[] CaculateFitness(List<List<List<Integer>>> firstPopulation,int VertexAllNumber,int Graphset[][][],
                                        int VertexPreNumber[],boolean VisitedNumbersSet[][])
    {
        int VertexNumberAfterFusion[]=new int[PopulationSize];
        Arrays.fill(VertexNumberAfterFusion,0);
        int EdgeNumberAfterFusion[]=new int[PopulationSize];
        int fitness[]=new int[PopulationSize];
        int populationIndex,row,col,graphIndex,i,j,x,y;

        //VisitedNumbersSet=caculateVistitedNumberSet(firstPopulation,VertexAllNumber,VertexPreNumber);
        //计算该种群每个编码生成点的数目
        VertexNumberAfterFusion=CaculateVertexNumberAfterFusion(firstPopulation,VertexAllNumber,
                Graphset,VertexPreNumber,VisitedNumbersSet);
        EdgeNumberAfterFusion=CaculateEdgeNumberAfterFusion(firstPopulation,VertexAllNumber,Graphset,
                         VertexPreNumber,VisitedNumbersSet);
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            fitness[populationIndex] = EdgeNumberAfterFusion[populationIndex];
                    //* VertexNumberAfterFusion[populationIndex];
        }

        return fitness;
    }

    //生成新个体
    public static List<List<Integer>> generateNewSingle(int Graphset[][][],int VertexAllNumber,int VertexPreNumber[])
    {
        int GraphIndex,VertexIndex,VertexNumberIndex;
        int firstPopulationSingleMatrix[][]=new int[GraphNumber][VertexAllNumber];
        for(GraphIndex=0;GraphIndex<GraphNumber;GraphIndex++)
        {
            for(VertexIndex=0;VertexIndex<VertexAllNumber;VertexIndex++)
            {
                firstPopulationSingleMatrix[GraphIndex][VertexIndex]=-1;
            }
        }

        Random random=new Random();
        for(GraphIndex=0;GraphIndex<GraphNumber;GraphIndex++)
        {
            //随机找到编码矩阵每一行中没有映射关系的位置,用当前节点值填充
            for(VertexNumberIndex=0;VertexNumberIndex<Graphset[GraphIndex].length;VertexNumberIndex++)
            {
                int pos=random.nextInt(VertexAllNumber);
                while(firstPopulationSingleMatrix[GraphIndex][pos]!=-1)
                {
                    pos=random.nextInt(VertexAllNumber);
                }
                firstPopulationSingleMatrix[GraphIndex][pos]=VertexNumberIndex;
            }

        }

        List<List<Integer>> firstpopulationSingle=new ArrayList<>();
        int row,col,index;
        //生成每个图(行)的对应数组,默认全为-1
        for(index=0;index<GraphNumber-1;index++)
        {
            List<Integer> graphMap=new ArrayList<>();
            for(int count=0;count<VertexPreNumber[index+1];count++)
            {
                graphMap.add(-1);
            }
            firstpopulationSingle.add(graphMap);
        }

        //将全矩阵编码映射成下三角
        for(col=0;col<firstPopulationSingleMatrix[0].length;col++)
        {
            int belong_col=-1;
            for(row=0;row<GraphNumber;row++)
            {
                //已经有了归属列
                if(belong_col>=0)
                {
                    //将全编码矩阵的row映射到半编码矩阵的row-1
                    //System.out.println("belong_col: "+belong_col);
                    //System.out.println("set Number: "+firstPopulationMatrix[populationIndex][row][col]);
                    firstpopulationSingle.get(row-1).set(belong_col,firstPopulationSingleMatrix[row][col]);
                }
                //如果没有归属列
                else if(belong_col==-1)
                {
                    if(firstPopulationSingleMatrix[row][col]==-1)
                    {
                        continue;
                    }
                    //如果不为-1,计算belong_col
                    else if(firstPopulationSingleMatrix[row][col]>-1)
                    {
                        belong_col=VertexPreNumber[row]+firstPopulationSingleMatrix[row][col];

                    }
                }

            }
        }
        return firstpopulationSingle;
    }


    //全矩阵编码转为下三角矩阵编码
    public static List<List<List<Integer>>> FullMatrixToLadderMatrix(int firstPopulationMatrix[][][],int VertexPreNumber[],int Graphset[][][],
                                                                     int VertexAllNumber)
    {
        List<List<List<Integer>>> firstpopulation=new ArrayList<>();
        int populationIndex=0;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            // System.out.println("population~~~~~~~~~~"+populationIndex);
            List<List<Integer>> firstpopulationSingle=generateFirstPopulationSingle(firstPopulationMatrix,VertexPreNumber,Graphset,
                    VertexAllNumber,populationIndex);
            firstpopulation.add(firstpopulationSingle);
        }
        return firstpopulation;
    }

    //生成下三角矩阵单体
    public static List<List<Integer>> generateFirstPopulationSingle(int firstPopulationMatrix[][][],int VertexPreNumber[],int Graphset[][][],
                                                                    int VertexAllNumber,int populationIndex)
    {
        List<List<Integer>> firstpopulationSingle=new ArrayList<>();
        int row,col,index;
        //生成每个图(行)的对应数组,默认全为-1
        for(index=0;index<GraphNumber-1;index++)
        {
            List<Integer> graphMap=new ArrayList<>();
            for(int count=0;count<VertexPreNumber[index+1];count++)
            {
                graphMap.add(-1);
            }
            firstpopulationSingle.add(graphMap);
        }

        //将全矩阵编码映射成下三角
        for(col=0;col<firstPopulationMatrix[populationIndex][0].length;col++)
        {
            int belong_col=-1;
            for(row=0;row<GraphNumber;row++)
            {
                //已经有了归属列
                if(belong_col>=0)
                {
                    //将全编码矩阵的row映射到半编码矩阵的row-1
                    //System.out.println("belong_col: "+belong_col);
                    //System.out.println("set Number: "+firstPopulationMatrix[populationIndex][row][col]);
                    firstpopulationSingle.get(row-1).set(belong_col,firstPopulationMatrix[populationIndex][row][col]);
                }
                //如果没有归属列
                else if(belong_col==-1)
                {
                    if(firstPopulationMatrix[populationIndex][row][col]==-1)
                    {
                        continue;
                    }
                    //如果不为-1,计算belong_col
                    else if(firstPopulationMatrix[populationIndex][row][col]>-1)
                    {
                        belong_col=VertexPreNumber[row]+firstPopulationMatrix[populationIndex][row][col];

                    }
                }

            }
        }

        return firstpopulationSingle;
    }


    public static int findColInGraphMatrix(int firstPopulationSingleMatrix[][],int number,int GraphIndex )
    {
        for(int index=0;index<firstPopulationSingleMatrix[GraphIndex].length;index++)
        {
            if(firstPopulationSingleMatrix[GraphIndex][index]==number)
            {
                return index;
            }
        }
        return 0;
    }

    //生成初始化种群,全矩阵编码,行数代表对应的图,列数为所以图节点之和
    public static int[][][] generateFirstPopulaitonMatrix(int Graphset[][][],int VertexAllNumber)
    {
        int firstPopulation[][][]=new int[PopulationSize][GraphNumber][VertexAllNumber];
        //对全矩阵初始化值全为-1,即都无映射关系
        int populationIndex,GraphIndex,VertexIndex,VertexNumberIndex;
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            for(GraphIndex=0;GraphIndex<GraphNumber;GraphIndex++)
            {
                for(VertexIndex=0;VertexIndex<VertexAllNumber;VertexIndex++)
                {
                    firstPopulation[populationIndex][GraphIndex][VertexIndex]=-1;
                }
            }
        }


        System.out.println("初始化赋值完成");
        //随机将每个图的每个节点映射到图编码矩阵中
        Random random=new Random();
        for(populationIndex=0;populationIndex<PopulationSize;populationIndex++)
        {
            for(GraphIndex=0;GraphIndex<GraphNumber;GraphIndex++)
            {
                //随机找到编码矩阵每一行中没有映射关系的位置,用当前节点值填充
                for(VertexNumberIndex=0;VertexNumberIndex<Graphset[GraphIndex].length;VertexNumberIndex++)
                {
                    int pos=random.nextInt(VertexAllNumber);
                    while(firstPopulation[populationIndex][GraphIndex][pos]!=-1)
                    {
                        pos=random.nextInt(VertexAllNumber);
                    }
                    firstPopulation[populationIndex][GraphIndex][pos]=VertexNumberIndex;
                }

            }
        }


        // System.out.println("firstPopulationParent "+Arrays.deepToString(firstPopulation[0]));
        //System.out.println("firstPopulationParent "+Arrays.deepToString(firstPopulation[25]));
        //System.out.println("firstPopulationParent "+Arrays.deepToString(firstPopulation[49]));
        return firstPopulation;
    }



    public static void GenerateGraph(String Filename,int Graphset[][][])
    {
        //逐行读取文件
        File file=new File(Filename);
        BufferedReader reader=null;
        try {
            System.out.println("以行为单位逐行读取");
            reader=new BufferedReader(new FileReader(file));
            int i;
            for(i=0;i<GraphNumber;i++)
            {
                //System.out.println(i);
                String line;
                int vertex=0;
                line=reader.readLine();//first line to show graph number
                //图顶点数
                while((line=reader.readLine()).charAt(0)=='v')
                {

                    String lineArray[]=line.split(" ");
                    vertex=Integer.parseInt(lineArray[1]);
                    //System.out.println(line+"line"+"vertex"+vertex);
                }
                //System.out.println(vertex+"vertex");
                Graphset[i]=new int[vertex+1][vertex+1];
                //Graphset[i][1][1]=2;

                if(line==null)//无边图
                {
                    continue;
                }
                else //有边图
                {
                    String lineArray[]=line.split(" ");
                    int vertex1=Integer.parseInt(lineArray[1]);
                    int vertex2=Integer.parseInt(lineArray[2]);
                    Graphset[i][vertex1][vertex2]=1;
                    Graphset[i][vertex2][vertex1]=1;
                    while(!(line=reader.readLine()).equals("endGraph"))
                    {

                        //System.out.println("line is: "+line);
                        lineArray=line.split(" ");
                        /*
                        if(line.charAt(0)=='t')
                        {
                            break;
                        }
                        */
                        if(lineArray.length>1)
                        {
                            vertex1=Integer.parseInt(lineArray[1]);
                            vertex2=Integer.parseInt(lineArray[2]);
                            //System.out.println("vertex1 "+vertex1);
                            //System.out.println("vertex2 "+vertex2);
                            Graphset[i][vertex1][vertex2]=1;
                            Graphset[i][vertex2][vertex1]=1;
                            //System.out.println("loop finish read edge");
                        }

                    }
                    //System.out.println("finish read edge");
                }
            }
            /*
            for(i=0;i<Graphset[9999].length;i++)
            {
                System.out.println(Arrays.deepToString(Graphset[9999]));
            }
            */
            System.out.println(Arrays.deepToString(Graphset[0]));
            reader.close();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }


    }
}
