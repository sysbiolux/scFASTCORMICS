# model_expansion.py
"""Python module demonstrates passing MATLAB types to Python functions"""
"""https://nl.mathworks.com/help/matlab/matlab_external/call-user-defined-custom-module.html"""
def search(words):
    """Return list of words containing 'son'"""
    newlist = [w for w in words if 'son' in w]
    return newlist

def theend(words):
    """Append 'The End' to list of words"""
    words.append('The End')
    return words


def display_model(model):
    """Append 'The End' to list of words"""
    print("Hello world")
    %print(model.summary())
    return