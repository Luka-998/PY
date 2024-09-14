import random
def print_box(text):
        print(len(text)*'*')
        print(text)
        print(len(text)*'*')
print_box("This is a game :)")


mm_count = random.randint(1,100)
attempt_limit = 3
attempts = 0

while attempts < attempt_limit:
    guess_number = input("How many M&M's are left?")
    attempts +=1
    guess = int(guess_number)
    if mm_count == guess_number:
        print(f"Good Job! You guessed correct! It is {guess_number}")
        break

    elif guess < mm_count:
        print("Sorry that is too low number")
    elif guess > mm_count:
        print ("Sorry that is too high")



