#include <stdio.h>

int main(int argc, char**argv) {

    // Initialize host variables

    // Allocate device variables

    // Copy tableau to Device

    // Iterate until optimal is found

        // Copy first line of the Simplex tableau to host

        // Find index of entering variable k

        // Copy index k to device

        // Kernel 1 to process ratio column

        // Copy ratio column to host

        // Find the index of the leaving variable r

        // Copy index r to device

        // Kernel 2 to update the line r of the Simplex tableau

        // Kernel 3 to update Simplex tableau

        // Kernel 4 to Update column k of the Simplex Tableau

        // Check if should stop

    // Calculate optimal value and return it

}
