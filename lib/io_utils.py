def set_test_mode(byc):

    test_mode = False

    try:
        if byc["args"].test:
            test_mode = True
            print( "¡¡¡ TEST MODE - no db update !!!")
    except:
        pass
        
    return test_mode