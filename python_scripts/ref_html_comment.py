with open('temp_lit_review.txt') as file:
    data = file.read()

first_bracket = '['
second_bracket = ']'

data.replace(first_bracket, '<!-- ' + first_bracket)

data.replace(second_bracket, second_bracket + ' -->')

with open('temp_lit_review.txt', 'w') as f:
    f.write(data)
