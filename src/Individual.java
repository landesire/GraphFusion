import java.util.ArrayList;
import java.util.List;

/**
 * Created by daibinding on 17/1/4.
 */
public class Individual {
    List<List<Integer>> graphFusionCode;
    int fitness;
    boolean VisitedNumberSetSingle[];

    //重载构造器
    Individual(List<List<Integer>> graphFusionCode,int fitness,boolean VisitedNumberSetSingle[])
    {
        this.graphFusionCode=new ArrayList<>();
        int index;
        //myAddALL方法,不涉及引用问题
        for( index=0;index<graphFusionCode.size();index++)
        {
            List<Integer> tempList=new ArrayList<>();
             tempList.addAll(graphFusionCode.get(index));
            this.graphFusionCode.add(tempList);
        }

        //this.clone()
        //this.graphFusionCode.addAll(graphFusionCode);
        this.fitness=fitness;
        this.VisitedNumberSetSingle=new boolean[VisitedNumberSetSingle.length];
        for(index=0;index<VisitedNumberSetSingle.length;index++ )
        {
            this.VisitedNumberSetSingle[index]=VisitedNumberSetSingle[index];
        }
    }


}
