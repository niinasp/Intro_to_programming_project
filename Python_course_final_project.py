#====================================================================================================================================================
# Sequences:
#====================================================================================================================================================
  
def askSeq(seq):
    """Checks a nucleotide sequence. Returns the sequence in lowercase letters if no errors are detected, otherwise returns True."""
 
    if len(seq) < 2:
        print ""
        print "The sequence is either missing or too short. Please enter a longer sequence.\n"
        return True
     
    # Check that the sequence only contain characters "a", "t", "c" or "g"
 
    nucleotides = 'atcgATCG'
     
    if len(seq) >= 2:
 
        # If non-nucleotide character is found, an error is printed out
        # Otherwise the input sequence modified to be in lowercese is returned
         
        for i in seq:
             
            if i not in nucleotides:
 
                print ""
                print "Non-nucleotide characters were detected in the sequence. Please enter a nucleotide sequence.\n"
                return True
 
            else:
                pass

    modSeq = str.lower(seq)
    return modSeq
 
  
  
#====================================================================================================================================================
# Parameters:
#====================================================================================================================================================
   
def defineSettings(dic):
    """Creates a dictionary of user defined substitution scores and also asks for gap penalty.
    Returns the dictionary and the gap penalty."""
       
    # Asks user to define settings to apply for match, mismatch and gap
   
    while True:
   
        try:
            match = input("Please enter a 'Match' score: ")
        except NameError:
            print "Please enter a number.\n"
        else:
            break
       
    while True:
        try:
            mismatch = input("Please enter a 'Mismatch' score: ")
        except NameError:
            print "Please enter a number.\n"
        else:
            break
   
    while True:
        try:
            userGp = input("Please enter a score for gap penalty: ")
        except NameError:
            print "Please enter a number.\n"
        else:
            break
       
    # Create a score dictionary based on user defined scores for match/mismatch
    # Use the same keys as in the dictionary that specifies the default settings (dictionary called "scores")
   
    userScores = {}
       
    for key in dic:
        if key == 'aa' or key == 'cc' or key == 'gg' or key == 'tt':
            userScores[key] = match  #for any keys where the two nucleotides are the same, the user defined match score is applied as the value
        else:
            userScores[key] = mismatch #in all other instances, the user defined mismatch score is applied as the value
               
    return userScores, userGp
   
    
#===================================================================================================================================================
#  Initialization Score Matrix
#===================================================================================================================================================
    
def initializeScoreMatrix(seq1, seq2, gp):
    """Initialises the score matrix for the input sequences using the chosen substitution scores and gap penalty.
    Returns the matrix."""
   
    # n: the number of rows in the matrix
    # m: the number of columns in the matrix
   
    n = len(seq2) + 2
    m = len(seq1) + 2
       
    # Create a matrix with n rows and m columns
    # Fill the matrix with hyphen characters
   
    matrix = [["-" for x in xrange(m)] for y in xrange(n)]
   
    # Fill the first row with characters of seq1
    # Start from the end of the row (i.e. the last column), iterate throuhg reversed seq1
       
    column = m - 1
       
    for i in seq1[::-1]:
        matrix[0][column] = i
        column = column - 1
   
    # Fill the first column with characters of seq2
    # Start from the bottom of the column (i.e. the last row), iterate throuhg reversed seq2
   
    row = n - 1
       
    for i in seq2[::-1]:
        matrix[row][0] = i
        row = row - 1
   
    # The starting position of the alignment (second row, second column) is filled with zero
   
    matrix[1][1] = 0
   
    # Fill the second row with gap penalties
    # Start from the third column in the second row and move towards left
    # When iterating over xrange(2, m)
    # the number 2 specifies the index of the starting column
    # and m refers to the number of columns in the matrix (the index of the last column is m - 1)
    # With each move, add one gap penalty to the score that is stored in the previous cell
   
    for i in xrange(2, m):
        matrix[1][i] = matrix[1][i - 1] + gp
           
    # For each cell, store the score and direction "L" as a tuple
   
    for i in xrange(2, m):
        matrix[1][i] = (matrix[1][i], "L")
           
    # Similarly, fill the second column with gap penalties in the matrix
    # Start from the third row in the second column and move downwards
    # When iterating over xrange(2, n)
    # the number 2 specifies the index of the starting row
    # and n refers to the number of rows in the matrix (the index of the last row is n - 1)
    # With each move, add one gap penalty to the score that is stored in the previous cell
   
    for j in xrange(2, n):
        matrix[j][1] = matrix[j - 1][1] + gp
   
    # For each cell, store the score and direction "U" as a tuple
   
    for j in xrange(2, n):
        matrix[j][1] = (matrix[j][1], "U")
   
    return matrix
    
      
#==================================================================================================================================================
# Global Pairwise Alignment Score Matrix
#==================================================================================================================================================
    
def calculateScoreMatrix(matrix, dic, gp):
    """Calculates the highest score for each cell in the matrix using the chosen substitution scores and gap penalty.
    Returns the matrix with the final scores."""
 
    # n: the number of rows in the score matrix
    # m: the number of columns in the score matrix
 
    n = len(matrix)
    m = len(matrix[0])
 
    # Iterate over the matrix row by row, column by column starting from the position [2][2]
 
    for i in xrange(2, n):
        for j in xrange(2, m):
             
            # For each cell, calculate three values (a, b, c) from existing scores in the upper, left and diagonal cells
 
            # When score is calculated from the upper cell ("U"), add gap penalty to the existing score
             
            a = matrix[i - 1][j][0] + gp
 
            # When score is calculated from the cell on the left ("L"), add gap penalty to the existing score
             
            b = matrix[i][j - 1][0] + gp
 
            # When score is calculated from the diagonal cell ("D"), add correct substitution score to the existing score
            # Substitutuion scores are read from the dictionary
            # Create the key for the dictionary by combining the nucleotides from the first row and first column that correspond to that cell
            # Get the substitution score from the dictionary and add that value to the existing score
 
            key = matrix[0][j] + matrix[i][0]
 
            # The position [2][2] does not contain a tuple, instead it contains 0
 
            if matrix[i - 1][j - 1] == 0:
                c = matrix[i - 1][j - 1] + dic[key]
 
            # Every other cell contains a tuple
            # Use index zero to get the score from the tuple
 
            else:
                c = matrix[i - 1][j - 1][0] + dic[key]
 
            # Choose the highest score
            # Store it together with the move ("U", "L", "D") made to reach the maximum score as a tuple
            # If there are two or three directions that yielded the highest score, store all directions to the tuple                
                 
            highest = max(a, b, c)
 
            if highest == a:
 
                if b == highest:
                    matrix[i][j] = (a, "U", "L")
                     
                elif c == highest:
                    matrix[i][j] = (a, "U", "D")
 
                elif b == highest and c == highest:
                    matrix[i][j] = (a, "U", "L", "D")
                     
                else:
                    matrix[i][j] = (a, "U")
 
            if highest == b:
                 
                if a == highest:
                    matrix[i][j] = (b, "U", "L")
                     
                if c == highest:
                    matrix[i][j] = (b, "L", "D")
             
                elif a == highest and c == highest:
                    matrix[i][j] = (b, "U", "L", "D")
                     
                else:
                    matrix[i][j] = (b, "L")
 
            if highest == c:
 
                if a == highest:
                    matrix[i][j] = (c, "U", "D")
 
                elif b == highest:
                    matrix[i][j] = (c, "L", "D")
 
                elif a == highest and b == highest:
                     matrix[i][j] = (c, "U", "L", "D")
                     
                else:
                    matrix[i][j] = (c, "D")
 
    return matrix
 
  
#===================================================================================================================================================
# Results
#===================================================================================================================================================
   
def outputResults(matrix, seq1, seq2, gp):
    """Outputs the alignment results to file named results.txt"""
      
    # Get score from the lowest right cell of score matrix
    score = matrix[len(matrix)-1][len(matrix[0])-1][0]
      
    # Get time of the analysis
    timestamp = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    # Open file for writing
    res = open("results.txt", "w")
 
  
    # Output results
    res.write("The results of the global pairwise alignment using dynamic programming \n \n")
    res.write("The analysis was run: "+ timestamp + "\n\n")
    res.write("--------------------------------------------" + "\n\n")
    res.write("Sequence 1: " + seq1 + "\n")
    res.write("Sequence 2: " + seq2 + "\n\n")
    res.write("Gap penalty used: " + str(gp) + "\n\n")
    res.write("--------------------------------------------" + "\n\n")
    res.write("Score of the best alignment found is: " + str(score) + "\n\n")
    res.write("--------------------------------------------" + "\n\n")
    res.write("Score matrix:" + "\n\n")
  
    # Use nested for loops to modify each cell for writing in file
    for i in matrix:
  
        # initialise empty matrix 
        e = []
        for j in i:
  
            # make length of each cell 15 (so that columns are aligned in the results matrix)
            j=str(j).ljust(15, " " )
  
            # append each modified cell into empty list
            e.append(j)
  
            # join list elements to string (one row in the matrix)
            line = "".join(e)+ "\n"
  
        # write row to file
        res.write(line)
    res.close()
  
    
    
#===================================================================================================================================================
# Main program
#===================================================================================================================================================
import datetime
print "Welcome to our nucleotide sequence alignment program"
print "This program creates a global pairwise alignmnet of user input sequences"
print ""
 
while True:
     
    # query for nucleotide sequences
   
    while True:
        sequence1 = raw_input("Please enter your first nucleotide sequence: ")
 
        seq1 = askSeq(sequence1)
        
        if isinstance(seq1, bool):
            continue
        else:
            break
 
    while True:
        sequence2 = raw_input("Please enter your second nucleotide sequence: ")
 
        seq2 = askSeq(sequence2)
     
        if isinstance(seq2, bool):
            continue
        else:
            break
   
    # notify the user of default scoring system and give the user option to select their own scoring system before proceeding
   
    print ""
    print "The Optimal global pairwise alignment score can be calculated one of two ways:"
    print ""  
    print "(1) Using default settings:  Match = 1, Mismatch = -1, Gap penalty = -1"
    print "(2) Using your own settings"
    print ""
  
    defScores = {'aa': 1, 'ac': -1, 'ag': -1, 'at': -1,
                 'ca': -1, 'cc': 1, 'cg': -1, 'ct': -1,
                 'ga': -1, 'gc': -1, 'gg': 1, 'gt': -1,
                 'ta': -1, 'tc': -1, 'tg': -1, 'tt': 1} 
    defGp = -1
  
  
    # Query for gap score & check for errors
  
    option1 = raw_input("Do you want to use default settings (Y/N)? ")
  
    # Series of acceptable answers for the default settings query
  
    yes = set(['yes', 'y', 'ye', 'es', 'y es', 'ye s'])
    no  = set(['no', 'n', 'n o', 'o'])
   
   
    # If "yes", use the default settings to initialise and calculate the final score matrix
  
    choice = None
    while not choice:
  
        choice = option1.lower()
  
        if choice in yes:
            a = initializeScoreMatrix(seq1, seq2, defGp)
            b = calculateScoreMatrix(a, defScores, defGp)
            print ""
   
    # If "no", ask for settings
    # Use these user-defined settings to initialise and calculate the final score matrix
  
        elif choice in no:
            print ""
            userScores, userGp = defineSettings(defScores)
            a = initializeScoreMatrix(seq1, seq2, userGp)
            b = calculateScoreMatrix(a, userScores, userGp)
            print ""
  
     
    #if no answer given or any other character used the following message will be posted
    #and the user will be asked again if they want to use default settings
        else:
            print ("Please respond with either a yes or no")
            option1 = raw_input("Do you want to use default settings (Y/N)?")
            choice = None
     
     
    # Select gap penalty (default or user input) for printing in the results file
    if choice == "n":
        gp = userGp
    else:
        gp = defGp
 
    # Make results file
    outputResults(b, seq1, seq2, gp)

    print "The results were saved in a text file results.txt."
 
    print ""
    answer = raw_input("Do you want to continue to make an another alignment (y/n)? \nType y to continue and any other character to quit: ")
    print ""
 
    if answer == "y":
        pass
    else:
        print "Bye!"
        break
        
