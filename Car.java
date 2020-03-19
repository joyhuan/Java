public abstract class Car {
	public abstract void run();
	
	public static void main(String args[]){
//		Car b= new Audi();    //upcasting
//		b.run();
		
		Car[] v = new Car[3];
		v[0] = new Audi(); 
		v[0].run();
		v[1] = new Audi(); 
		v[1].run();
		v[0].run();
	}
}

class Audi extends Car {
	public void run(){
		System.out.println("Audi is running safely with 100km");
	}
	private int speed = 100; 
	public int getSpeed() {
		return speed; 
	}
}

