
public class CheckedDemo {
	/** Method main can relieve itself of this responsibility by having 
	 * its own throws clause
	 * The Java runtime system now has the responsibility of catching 
	 * Exceptions, since it calls main, and it will do so. RuntimeExceptions
	 * and Errors are not checked. 
	 * @param pars
	 * @throws Exception
	 */
	public static void main(String[] pars) throws Exception {
		first(); 
	}
	
	/** The occurrence of such a throws clause in the header
	of a method definition relieves the method of the responsibility 
	of catching objects of the mentioned classes and places that 
	burden on any method that calls this one. 
	*/ 
	public static void first() throws Exception{
		throw new Exception(); 
	}
}
