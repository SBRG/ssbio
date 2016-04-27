import datetime

class Date():
    def __init__(self):
        self.short_date = self.date_prefix()

    def date_prefix(self):
        today = datetime.date.today()
        return today.strftime('%y%m%d')

if __name__ == '__main__':
    d = Date()
    print(d.short_date)