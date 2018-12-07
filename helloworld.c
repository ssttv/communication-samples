#include <stdio.h>
int main()
{
    int sumOfTwoNumbers;
    int firstNumber;
    int secondNumber;
    
    printf("Enter two integers: \n");
 
    // Two integers entered by user is stored using scanf() function
    scanf("%d %d", &firstNumber, &secondNumber);
 
    // sum of two numbers in stored in variable sumOfTwoNumbers
    sumOfTwoNumbers = firstNumber + secondNumber;
 
    // Displays sum      
    printf("%d + %d = %d \n", firstNumber, secondNumber, sumOfTwoNumbers);
 
    return 0;
}