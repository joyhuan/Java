import static org.junit.Assert.*;

import org.junit.Test;

/** Contains methods to test all the methods in class ComparableArrays.
@author David Gries */
public class ComparableArraysTest {

    @Test
    public void testLinearSearch() {
        Integer[] b= {1, 3, 3, 3, 3, 4, 4, 6, 7, 8, 8, 10};
        assertEquals(12, ComparableArrays.linearSearch(b, -5));
        assertEquals(0, ComparableArrays.linearSearch(b, 1));
        assertEquals(12, ComparableArrays.linearSearch(b, 2));
        assertEquals(9, ComparableArrays.linearSearch(b, 8));
        assertEquals(12, ComparableArrays.linearSearch(b, 12));
    }

    @Test
    public void testBinarySearch() {
        Integer[] b= {1, 3, 3, 3, 3, 4, 4, 6, 7, 8, 8, 10};
        assertEquals(-1, ComparableArrays.binarySearch(b, -5));
        assertEquals(0, ComparableArrays.binarySearch(b, 1));
        assertEquals(0, ComparableArrays.binarySearch(b, 2));
        assertEquals(4, ComparableArrays.binarySearch(b, 3));
        assertEquals(6, ComparableArrays.binarySearch(b, 4));
        assertEquals(6, ComparableArrays.binarySearch(b, 5));
        assertEquals(10, ComparableArrays.binarySearch(b, 9));
        assertEquals(11, ComparableArrays.binarySearch(b, 10));
        assertEquals(11, ComparableArrays.binarySearch(b, 11));
    }

    @Test
    public void testMin() {
        Integer[] b= {3, 8, 5, -2};
        assertEquals(3, ComparableArrays.min(b, 0, 3));
        assertEquals(1, ComparableArrays.min(b, 1, 1));
        assertEquals(0, ComparableArrays.min(b, 0, 2));
    }

    @Test
    public void testSelectionSort() {
        Integer[] b= {};
        ComparableArrays.selectionSort(b);
        assertEquals("[]", ComparableArrays.toString(b));

        b= new Integer[] {6};
        ComparableArrays.selectionSort(b);
        assertEquals("[6]", ComparableArrays.toString(b));

        b= new Integer[] {6, 3, -2, 7, 5, 8, 7};
        ComparableArrays.selectionSort(b);
        assertEquals("[-2, 3, 5, 6, 7, 7, 8]", ComparableArrays.toString(b));
    }

    @Test
    public void testSelectionSort1() {
        Integer[] b= {};
        ComparableArrays.selectionSort1(b);
        assertEquals("[]", ComparableArrays.toString(b));

        b= new Integer[] {6};
        ComparableArrays.selectionSort1(b);
        assertEquals("[6]", ComparableArrays.toString(b));

        b= new Integer[] {6, 3, -2, 7, 5, 8, 7};
        ComparableArrays.selectionSort1(b);
        assertEquals("[-2, 3, 5, 6, 7, 7, 8]", ComparableArrays.toString(b));
    }

    @Test
    public void testInsertValue() {
        Integer[] b= {2, 4, 6, 7, 8, 5, 1};
        ComparableArrays.insertValue(b, 1, 5);
        assertEquals("[2, 4, 5, 6, 7, 8, 1]", ComparableArrays.toString(b));

        ComparableArrays.insertValue(b, 1, 6);
        assertEquals("[2, 1, 4, 5, 6, 7, 8]", ComparableArrays.toString(b));

        ComparableArrays.insertValue(b, 0, 1);
        assertEquals("[1, 2, 4, 5, 6, 7, 8]", ComparableArrays.toString(b));
    }

    @Test
    public void testInsertionSort() {
        Integer[] b= {};
        ComparableArrays.insertionSort(b, 0, -1);
        assertEquals("[]", ComparableArrays.toString(b));

        b= new Integer[] {6};
        ComparableArrays.insertionSort(b, 0, 0);
        assertEquals("[6]", ComparableArrays.toString(b));

        b= new Integer[] {6, 3, -2, 7, 5, 8, 7};
        ComparableArrays.insertionSort(b, 1, 6);
        assertEquals("[6, -2, 3, 5, 7, 7, 8]", ComparableArrays.toString(b));
    }

    @Test
    public void testPartition0() {
        Integer[] b= {4, 3, 8, 7, 8, 5, 1};
        int j= ComparableArrays.partition0(b, 0, 6);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 4);
        assertTrue(b[j] == 4);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 4);

        b= new Integer[]{0, 3, 8, 7, 8, 5, 1};
        j= ComparableArrays.partition0(b, 0, 6);
        assertEquals(0, j);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 0);
        assertTrue(b[j] == 0);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 0);
    }

    @Test
    public void testPartition() {
        Integer[] b= {4, 3, 8, 7, 8, 5, 1};
        int j= ComparableArrays.partition(b, 0, 6);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 4);
        assertTrue(b[j] == 4);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 4);

        b= new Integer[]{0, 3, 8, 7, 8, 5, 1};
        j= ComparableArrays.partition(b, 0, 6);
        assertEquals(0, j);
        for (int k= 0; k < j; k= k+1)
            assertTrue(b[k] < 0);
        assertTrue(b[j] == 0);
        for (int k= j+1; k < b.length; k= k+1)
            assertTrue(b[k] > 0);
    }

    @Test
    public void testMedianOf3() {
        Integer[] b= {4, 3, 8, 7, 8, 5, 1};
        ComparableArrays.medianOf3(b, 0, 6);
        assertEquals("[4, 3, 8, 7, 8, 5, 1]", ComparableArrays.toString(b));

        b= new Integer[]{1, 3, 8, 7, 8, 5, 4};
        ComparableArrays.medianOf3(b, 0, 6);
        assertEquals("[4, 3, 8, 7, 8, 5, 1]", ComparableArrays.toString(b));

        b= new Integer[]{1, 3, 8, 4, 8, 5, 7};
        ComparableArrays.medianOf3(b, 0, 6);
        assertEquals("[4, 3, 8, 1, 8, 5, 7]", ComparableArrays.toString(b));
    }

    @Test
    public void testCopy() {
        Integer[] b= {4, 3, 8, 7, 8, 5, 1};
        Comparable[] c= ComparableArrays.copy(b, 2, 4);
        assertEquals("[8, 7, 8]", ComparableArrays.toString(c));

        c= ComparableArrays.copy(b, 2, 1);
        assertEquals("[]", ComparableArrays.toString(c));
    }

    @Test
    public void testMerge() {
        Integer[] b= {2, 4, 5, 7, 9, 10, 0, 2, 3, 4, 8};
        ComparableArrays.merge(b, 0, 5, 10);
        assertEquals("[0, 2, 2, 3, 4, 4, 5, 7, 8, 9, 10]", ComparableArrays.toString(b));

        b= new Integer[]{2, 4, 5};
        ComparableArrays.merge(b, 0, 2, 2);
        assertEquals("[2, 4, 5]", ComparableArrays.toString(b));
    }

    @Test
    public void testMergeSort() {
        Integer[] b= new Integer[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= b.length - k;
        }
        ComparableArrays.mergeSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertTrue(k+1 == b[k]);
        }

        b= new Integer[500];
        for(int k= 0; k < b.length; k= k+1) {
            b[k]= k;
        }
        ComparableArrays.mergeSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertTrue(k == b[k]);
        }
    }

    @Test
    public void testQuicksort0() {
        // Test on an array {500, 499, 498, 497, ..., 1}
        Integer[] b= new Integer[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= b.length - k;
        }
        ComparableArrays.quickSort0(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertTrue(k+1 == b[k]);
        }

        // Test on an array {0, 1, 2, ..., 499}
        b= new Integer[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= k;
        }
        ComparableArrays.quickSort0(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertTrue(k == b[k]);
        }
    }

    @Test
    public void testQuicksort() {
        // Test on an array {500, 499, 498, 497, ..., 1}
        Integer[] b= new Integer[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= b.length - k;
        }
        ComparableArrays.quickSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertTrue(k+1 == b[k]);
        }

        // Test on an array {0, 1, 2, ..., 499}
        b= new Integer[500];
        for (int k= 0; k < b.length; k= k+1) {
            b[k]= k;
        }
        ComparableArrays.quickSort(b, 0, b.length-1);
        for (int k= 0; k < b.length; k= k+1) {
            assertTrue(k == b[k]);
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
