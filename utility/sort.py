
def natural( l ): 
    """ 
    Sort the given list in the way that humans expect. For example,
    
    >>> l = ["1", "4", "2" , "0"]
    >>> sort_nicely(l)
    >>> l
    ['0', '1', '2', '4']
    
    Note that negative values (its argument takes only string!) is shifted to the backward  

    >>> l = ["-1","3","-2", "4"]
    >>> sort_nicely(l)
    >>> l
    ['3', '4', '-1', '-2']
     
    """ 
    import re     
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key )    