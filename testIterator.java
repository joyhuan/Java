import java.util.ArrayList;
import java.util.Iterator;


public class testIterator {
	public static void main(String[] args) {
		ArrayList<String> list = new ArrayList<String>(); 
		list.add("A");
		list.add("B");
		list.add("C");
		list.add("D");
		
		Iterator iterator = list.iterator();

		System.out.println("List elements : ");
		
		while(iterator.hasNext())
			System.out.print(iterator.next() + " ");
		System.out.println();
	}
}
