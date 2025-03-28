#password generator 1) number of chars
#2) lenght

import random

characters = "sadpkoj231eo123sadsdik1-20ekaspocasz//12!"
number = input('Number of passwords: ')
number = int(number)
lenght = input("Input your password lenght: ")
lenght = int(lenght)
print ("\nHere are your passwords: ")
for pwd in range(number):
    passwords = ''
    for c in range(lenght):
        passwords += random.choice(characters)

    print(passwords)