def yes_no(question):
   while True:
      x = input(question + "Answer Yes/No only.").strip()
      if x == "Yes":
         return True
      elif x == "No":
         return False
      else:
         print("Please only type Yes/No.")

if yes_no("Do you like ice cream?"):
    print("like")
else:
    print("don't like")