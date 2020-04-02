#include <iostream>
int main()
{
  int input_var = 0;
  // Enter the do while loop and stay there until either
  // a non-numeric is entered, or -1 is entered.  Note that
  // cin will accept any integer, 4, 40, 400, etc.
  do {
    std::cout << "Enter a number (-1 = quit): ";
    // The following line accepts input from the keyboard into
    // variable input_var.
    // cin returns false if an input operation fails, that is, if
    // something other than an int (the type of input_var) is entered.
    if (!(cin >> input_var)) {
      std::cout << "Please enter numbers only." << std::endl;
      cin.clear();
      cin.ignore(10000,'\n');
    }
    if (input_var != -1) {
      std::cout << "You entered " << input_var << std::endl;
    }
  }
  while (input_var != -1);
  std::cout << "All done." << std::endl;

  return 0;
}