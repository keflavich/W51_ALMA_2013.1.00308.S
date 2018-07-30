__rethrow_casa_exceptions = True
h_init()
try:
    hifa_restoredata (vis=['uid___A002_Xc56496_X1912', 'uid___A002_Xc5a30f_X8c3', 'uid___A002_Xc5b7d7_X870', 'uid___A002_Xc5b7d7_X492d', 'uid___A002_Xc5b7d7_X4cdb'], session=['session_1', 'session_4', 'session_6', 'session_7', 'session_7'], ocorr_mode='ca')
finally:
    h_save()
