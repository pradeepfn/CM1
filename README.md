Modified CM1 code to work with phoenix checkpoint restart lib.


run
----



restart
--------


cm1.F contains the called to restart_read. We have to allocate the heap memory for
those variables using our checkpoint lib.
