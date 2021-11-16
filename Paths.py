import datetime
import os


def get_formated_date(date: datetime.datetime, format: str) -> str:
    """
    Returns the date in a specified format
    :param date: date
    :param format: format string
    :return: formated date
    """
    return date.strftime(format)


def create_folder(path: str) -> None:
    """
    Creates a new folder
    :param path: path
    :return: None
    """
    try:
        os.makedirs(path)
    except:
        pass


class Directories:

    def __init__(self, type: str):
        self.today = datetime.datetime.now()
        self.code_dir = os.path.dirname(os.path.abspath(__file__))
        self.output_folder = os.path.join(self.code_dir, "Results\\" + type + "\\"+ get_formated_date(self.today, "%Y") + "\\" \
                                          + get_formated_date(self.today, "%m") + "\\" +
                                          get_formated_date(self.today, "%Y_%m_%d") + "\\Time_" +
                                          get_formated_date(self.today, "%H_%M") + "\\")