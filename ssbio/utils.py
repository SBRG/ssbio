import datetime

def force_string(val=None):
    """Force a string representation of an object

    Args:
        val: object to parse into a string

    Returns:

    """
    if val is None:
        return ''
    return val if isinstance(val, str) else ';'.join(val)

def listn(val=None):
    """Force a list representation of an object

    Args:
        val: object to parse into a list

    Returns:

    """
    if val is None:
        return []
    return val if isinstance(val, list) else [val]

class Date():
    def __init__(self):
        self.short_date = self.date_prefix()

    def date_prefix(self):
        today = datetime.date.today()
        return today.strftime('%y%m%d')

if __name__ == '__main__':
    d = Date()
    print(d.short_date)