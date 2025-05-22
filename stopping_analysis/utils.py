def append_to_dict(dictionary, key, value):
    """
    appends to a dictionary, where the value is a list of elements

    Paramters:
    dictionary (dict): dictionary
    key (string): key
    value: (anything??): value to append to the list
    """

    if key not in dictionary:
        dictionary[key] = [value]
    else:
        dictionary[key].append(value)

    return dictionary
