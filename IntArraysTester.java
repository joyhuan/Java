import static org.junit.Assert.*;

import org.junit.Test;

/** Contains methods to test all the methods in class IntArrays.
@author David Gries */
public class IntArraysTester {
    
    @Test
    public void testLinearSearch() {
        int[] b= {1, 3, 3, 3, 3, 4, 4, 6, 7, 8, 8, 10};
        assertEquals(12, IntArrays.linearSearch(b, -5));
        assertEquals(0, IntArrays.linearSearch(b, 1));
        assertEquals(12, IntArrays.linearSearch(b, 2));
        assertEquals(9, IntArrays.linearSearch(b, 8));
        assertEquals(12, IntArrays.linearSearch(b, 12));
    }
    
    @Test
    public void testBinarySearch() {
        int[] b= {1, 3, 3, 3, 3, 4, 4, 6, 7, 8, 8, 10};
        assertEquals(-1, IntArrays.binarySearch(b, -5));
        assertEquals(0, IntArrays.binarySearch(b, 1));
        assertEquals(0, IntArrays.binarySearch(b, 2));
        assertEquals(4, IntArrays.binarySearch(b, 3));
        assertEquals(6, IntArrays.binarySearch(b, 4));
        assertEquals(6, IntArrays.binarySearch(b, 5));
        assertEquals(10, IntArrays.binarySearch(b, 9));
        assertEquals(11, IntArrays.binarySearch(b, 10));
        assertEquals(11, IntArrays.binarySearch(b, 11));
    }

    @Test
    public void testMin() {
        int[] b= {3, 8, 5, -2};
        assertEquals(3, IntArrays.min(b, 0, 3));
        assertEquals(1, IntArrays.min(b, 1, 1));
        assertEquals(0, IntArrays.min(b, 0, 2));
    }

    @Test
    public void testSelectionSort() {
        int[] b= {};
        IntArrays.selectionSort(b);
        assertEquals("[]", IntArrays.toString(b));

        b= new int[] {6};
        IntArrays.selectionSort(b);
        assertEquals("[6]", IntArrays.toString(b));

        b= new int[] {6, 3, -2, 7, 5, 8, 7};
        IntArrays.selectionSort(b);
        assertEquals("[-2, 3, 5, 6, 7, 7, 8]", IntArrays.toString(b));
    }
    
    @Test
    public void testSelectionSort1() {
        int[] b= {};
        IntArrays.selectionSort1(b);
        assertEquals("[]", IntArrays.toString(b));

        b= new int[] {6};
        IntArrays.selectionSort1(b);
        assertEquals("[6]", IntArrays.toString(b));

        b= new int[] {6, 3, -2, 7, 5, 8, 7};
        IntArrays.selectionSort1(b);
        assertEquals("[-2, 3, 5, 6, 7, 7, 8]", IntArrays.toString(b));
    }

    @Test
    public void testInsertValue() {
        int[] b= {2, 4, 6, 7, 8, 5, 1};
        IntArrays.insertValue(b, 1, 5);
        assertEquals("[2, 4, 5, 6, 7, 8, 1]", IntArrays.toString(b));
        
        IntArrays.insertValue(b, 1, 6);
        assertEquals("[2, 1, 4, 5, 6, 7, 8]", IntArrays.toString(b));
        
        IntArrays.insertValue(b, 0, 1);
        assertEquals("[1, 2, 4, 5, 6, 7, 8]", IntArrays.toString(b));
    }
    
    @Test
    public void testInsertionSort() {
        int[] b= {};
        IntArrays.insertionSort(b, 0, -1);
        assertEquals("[]", IntArrays.toString(b));

        b= new int[] {6};
        IntArrays.insertionSort(b, 0, 0);
        assertEquals("[6]", IntArrays.toString(b));

        b= new int[] {6, 3, -2, 7, 5, 8, 7};
        IntArrays.insertionSort(b, 1, 6);
        assertEquals("[6, -2, 3, 5, 7, 7, 8]", IntArrays.toString(b));
    }
    
    @Test
    public void testPartition0() {
        int[] b= {4, 3, 8, 7, 8, 5, 1};
        int j= IntArrays.partition0(b, 0, 6);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 4);
        assertTrue(b[j] == 4);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 4);
        
        b= new int[]{0, 3, 8, 7, 8, 5, 1};
        j= IntArrays.partition0(b, 0, 6);
        assertEquals(0, j);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 0);
        assertTrue(b[j] == 0);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 0);
    }
    
    @Test
    public void testPartition() {
        int[] b= {4, 3, 8, 7, 8, 5, 1};
        int j= IntArrays.partition(b, 0, 6);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 4);
        assertTrue(b[j] == 4);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 4);
        
        b= new int[]{0, 3, 8, 7, 8, 5, 1};
        j= IntArrays.partition(b, 0, 6);
        assertEquals(0, j);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 0);
        assertTrue(b[j] == 0);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 0);
    }
    
    @Test
    public void testMedianOf3() {
        int[] b= {4, 3, 8, 7, 8, 5, 1};
        IntArrays.medianOf3(b, 0, 6);
        assertEquals("[4, 3, 8, 7, 8, 5, 1]", IntArrays.toString(b));

        b= new int[]{1, 3, 8, 7, 8, 5, 4};
        IntArrays.medianOf3(b, 0, 6);
        assertEquals("[4, 3, 8, 7, 8, 5, 1]", IntArrays.toString(b));

        b= new int[]{1, 3, 8, 4, 8, 5, 7};
        IntArrays.medianOf3(b, 0, 6);
        assertEquals("[4, 3, 8, 1, 8, 5, 7]", IntArrays.toString(b));
    }
    
    @Test
    public void testCopy() {
        int[] b= {4, 3, 8, 7, 8, 5, 1};
        int[] c= IntArrays.copy(b, 2, 4);
        assertEquals("[8, 7, 8]", IntArrays.toString(c));
        
        c= IntArrays.copy(b, 2, 1);
        assertEquals("[]", IntArrays.toString(c));
    }
    
    @Test
    public void testMerge() {
        int[] b= {2, 4, 5, 7, 9, 10, 0, 2, 3, 4, 8};
        IntArrays.merge(b, 0, 5, 10);
        assertEquals("[0, 2, 2, 3, 4, 4, 5, 7, 8, 9, 10]", IntArrays.toString(b));
        
        b= new int[]{2, 4, 5};
        IntArrays.merge(b, 0, 2, 2);
        assertEquals("[2, 4, 5]", IntArrays.toString(b));
    }
    
    @Test
    public void testMergeSort() {
        int[] b= new int[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= b.length - k;
        }
        IntArrays.mergeSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertEquals(k+1, b[k]);
        }
        
        b= new int[500];
        for(int k= 0; k < b.length; k= k+1) {
            b[k]= k;
        }
        IntArrays.mergeSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertEquals(k, b[k]);
        }
    }
    
    @Test
    public void testQuickSort0() {
        // Test on an array {500, 499, 498, 497, ..., 1}
        int[] b= new int[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= b.length - k;
        }
        IntArrays.quickSort0(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertEquals(k+1, b[k]);
        }
        
     // Test on an array {0, 1, 2, ..., 499}
        b= new int[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= k;
        }
        IntArrays.quickSort0(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertEquals(k, b[k]);
        }
    }
    
    @Test
    public void testQuickSort() {
     // Test on an array {500, 499, 498, 497, ..., 1}
        int[] b= new int[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= b.length - k;
        }
        IntArrays.quickSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertEquals(k+1, b[k]);
        }
        
     // Test on an array {0, 1, 2, ..., 499}
        b= new int[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= k;
        }
        IntArrays.quickSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertEquals(k, b[k]);
        }
    }
    
    @Test
    public void testDutchNationalFlag() {
        // To make testing easier, take advantge of the fact that
        // it only distinguishes between <0, 0, and >0 and have
        //only 3 different values in the array
        int[] b= new int[] {-3, 5, 0, 0, 5, -3, 5, -3, -3, -3};
        IntArrays.DutchNationalFlag(b);
        assertEquals("[-3, -3, -3, -3, -3, 0, 0, 5, 5, 5]", IntArrays.toString(b));
        
        b= new int[] {5, 0, 0, 5, 5};
        IntArrays.DutchNationalFlag(b);
        assertEquals("[0, 0, 5, 5, 5]", IntArrays.toString(b));
        
        b= new int[] {5, 5, 5};
        IntArrays.DutchNationalFlag(b);
        assertEquals("[5, 5, 5]", IntArrays.toString(b));
        
        b= new int[] {};
        IntArrays.DutchNationalFlag(b);
        assertEquals("[]", IntArrays.toString(b));
        
        
    }
}
