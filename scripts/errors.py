
class OutOfSizeError(Exception):
    '''
    Error raised when one is trying to reach the data out of the size of file.
    '''

    def __init__(self, pointer):
        self.message = "Pointer on " + str(pointer) + " exceeds the size of data."
        super(Exception, self).__init__()
