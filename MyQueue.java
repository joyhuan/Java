import java.util.ArrayList;
import java.util.Stack;

public class MyQueue {
	Stack<Integer> que;
	
	/** Initialize your data structure here. */
	Stack<Integer> queue = new Stack<Integer>();
	// Push element x to the back of queue.
	
	// Only change is here 换杯子倒水的感觉
	public void push(int x) {
	    Stack<Integer> temp = new Stack<Integer>();
	    while(!queue.empty()){
	        temp.push(queue.pop());
	    }
	    queue.push(x);
	    while(!temp.empty()){
	        queue.push(temp.pop());
	    }
	}

	// Removes the element from in front of queue.
	public void pop() {
	    queue.pop();
	}

	// Get the front element.
	public int peek() {
	    return queue.peek();
	}

	// Return whether the queue is empty.
	public boolean empty() {
	    return queue.empty();
	}
}

/**
 * Your MyQueue object will be instantiated and called as such:
 * MyQueue obj = new MyQueue();
 * obj.push(x);
 * int param_2 = obj.pop();
 * int param_3 = obj.peek();
 * boolean param_4 = obj.empty();
 */