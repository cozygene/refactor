
import logging
from logging import handlers
import sys
import time
import traceback

LOGGING_LEVEL = logging.DEBUG
LOG_FILE_PATH = 'refactor.log'


def format_exc_info(exc_info):
    exceptionType, exception, tracebackObject = exc_info
    tracebackString = ' ---> '.join(['%s:%d, in %s: %s' % (filename, lineNumber, functionName, text) for (filename, lineNumber, functionName, text) in traceback.extract_tb(tracebackObject)])
    return '; EXCEPTION: %s, %s; TRACEBACK: %s' % (exceptionType, exception, tracebackString)

class _Formatter(logging.Formatter):
    def __init__(self, namespace, *args, **kwargs):
        self._namespace = namespace
        super(_Formatter, self).__init__(*args, **kwargs)

    def formatException(self, record):
        if record.exc_info is not None:
            return format_exc_info(record.exc_info)
        else:
            return ''

class _FileFormatter(_Formatter):
    def format(self, record):
        timestamp = self.formatTime(record)
        source = '{}:{}'.format(record.filename, record.lineno)
        message = record.getMessage()
        exception = self.formatException(record)
        
        return '{} glint {:<10} {:<10} {:<21} {}{}'.format(timestamp, self._namespace, record.levelname, source, message, exception)


class _ConsoleFormatter(_Formatter):
    def format(self, record):
        source = '{}:{}'.format(record.filename, record.lineno)
        message = record.getMessage()
        exception = self.formatException(record)
        
        return '{:<10} {}{}'.format(record.levelname, message, exception)



def configureLogging(namespace, syslog_host='localhost'):
    logging.raiseExceptions = 0
    logging.captureWarnings(True)

    logger = logging.getLogger()
    logger.setLevel(LOGGING_LEVEL)

    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(_ConsoleFormatter(namespace))
    streamHandler.setLevel(LOGGING_LEVEL)
    logger.addHandler(streamHandler)


    fileHandler = logging.FileHandler(LOG_FILE_PATH,mode='w')
    fileHandler.setLevel(LOGGING_LEVEL)
    fileHandler.setFormatter(_FileFormatter(namespace))
    logger.addHandler(fileHandler)