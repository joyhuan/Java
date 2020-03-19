
public class Person {
	private String firstName; //not null
	private String middleName;
	private String lastName;
	public Person(String f, String l) {
		assert f != null;
		firstName= f;
		lastName= l;
	}
	public Person(String f, String m, String l) {
		this(f, l);
		middleName= m;
	}
}
