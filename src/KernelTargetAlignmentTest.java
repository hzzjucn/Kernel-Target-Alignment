import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class KernelTargetAlignmentTest {

	public static void main(String[] args) {
		KernelTargetAlignment kta = new KernelTargetAlignment();
				
		RealMatrix t1 = new Array2DRowRealMatrix(new double[][] { 
				{ 1.68, 0.34, 0.38 },
				{ 0.34, 3.09, -1.59 }, 
				{ 0.38, -1.59, 1.54 } });
		
		RealMatrix t2 = new Array2DRowRealMatrix(new double[][] { 
				{ 1, 2, 3},
				{ 3, 3, -1}, 
				{ 4,-1, 8}});
		// test Center
//		System.out.println((kta.Center (t1)).toString());
		
		//test frobenius product
//		System.out.println(kta.FrobeniusProduct(t1, t2));
		
		//test normalization
		System.out.println (kta.Normalize(t2));
	}
}
