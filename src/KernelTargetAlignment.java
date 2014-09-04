import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class KernelTargetAlignment {		
	
	//get centered kernel matrix
	public RealMatrix Center (RealMatrix km)
	{
		double dim = km.getRowDimension();
		RealMatrix I = MatrixUtils.createRealIdentityMatrix((int)dim);
		System.out.println(I.toString());
		
		RealMatrix one = MatrixUtils.createRealMatrix((int)dim, (int)dim);
		
		for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			{
			   one.setEntry(i, j, 1d);
			}		
		System.out.println(one.toString());
		
		RealMatrix u = I.subtract(one.scalarMultiply(1d/dim));		
		System.out.println(u.toString());
		
	    return (u.multiply(km)).multiply(u);
	}
	
	//normalize kernels by: K(x,y)/sqrt(K(x,x)K(y,y)) 
	public RealMatrix Normalize (RealMatrix km)
	{
		int dim = km.getRowDimension();
		double[]diag = new double[dim];
		
		for (int k=0; k<dim; k++)
		{
		    diag[k] = km.getEntry(k, k);
		}		
		RealMatrix d = MatrixUtils.createColumnRealMatrix (diag);		
		System.out.println (d.multiply(d.transpose()));
		
		double[][] tmp = d.multiply(d.transpose()).getData();	
		for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			{
				tmp[i][j]= km.getEntry(i, j)/Math.sqrt(tmp[i][j]);
			}				
		return MatrixUtils.createRealMatrix(tmp);
	}
	
	// the frobenius product of two kernel matrixes
	public double FrobeniusProduct(RealMatrix km1, RealMatrix km2)
	{		
	   return (km1.transpose().multiply(km2)).getTrace();		
	}
		
	public static boolean checkIsPSD (RealMatrix km)
	{
		EigenDecomposition ed = new EigenDecomposition (km);
		double[] egValues = ed.getRealEigenvalues();
		
		double min = Double.MAX_VALUE;
		for (double e_val: egValues)
		{
			if (e_val < min)
				min = e_val;
		}
		if (min < -1.0E-5)
			return false;
		else
			return true;					
	}
		
    public static Object readMatrix (String path) throws IOException, ClassNotFoundException{         
        FileInputStream in = new FileInputStream(path);  
        ObjectInputStream objread = new ObjectInputStream(in); 
        Object map = objread.readObject();       
        objread.close(); 
        return map;
    }
    	
	public static void main(String[] args) throws ClassNotFoundException, IOException 
	{
		// read matrixes
		String folder = "data/";
		String[] modalities = {"name", "category", "dev", "description", "update", "permission", "image", "content", "size", "review"};
		int p = modalities.length;
		String target = "Y";		
		ArrayList <RealMatrix> kms = new ArrayList <RealMatrix> ();
		
		for (String modality: modalities)
		{
			RealMatrix km = (RealMatrix) readMatrix (folder + modality + ".data");
			kms.add (km);
		}		
		RealMatrix ky = (RealMatrix) readMatrix (folder + target + ".data");
		
		//Construct M, a
		RealMatrix a = MatrixUtils.createRealMatrix (p, 1);
		RealMatrix M = MatrixUtils.createRealMatrix (p, p);
		
		KernelTargetAlignment kta = new KernelTargetAlignment();
		
	 	RealMatrix ky_c = kta.Center(ky);
	 	ky_c = kta.Normalize(ky_c);
	 	
	 	ArrayList <RealMatrix> km_cs = new ArrayList <RealMatrix> ();	 	
	 	
	 	for (int i=0; i<kms.size(); i++)
	 	{
	 		RealMatrix km_c = kta.Center (kms.get(i));
	 		km_c = kta.Normalize (km_c);	 		
	 		km_cs.set (i, km_c);
	 		a.setEntry (i, 1, kta.FrobeniusProduct (km_c, ky_c));
	 	}
	 	
	 	for (int i=0; i<p; i++)
	 		for (int j=0; j<p; j++)
	 		{
	 			double entry = kta.FrobeniusProduct(km_cs.get(i), km_cs.get(j));
	 			M.setEntry(i, j, entry);
	 			M.setEntry(j, i, entry);	 				 			
	 		}
	 	
	 	// solve the QP optimization
	 	
	}
}
