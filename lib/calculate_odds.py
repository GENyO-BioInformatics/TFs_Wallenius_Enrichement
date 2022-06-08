




def getGAModds(lens,DEs):
    x = numpy.array(lens)
    y = numpy.array(DEs)
    ww = numpy.argsort(x) # ww
    size = math.ceil(len(y) / 10)
    low = sum(y[ww][0:size])
    hi = sum(y[ww][(len(y) - size):len(y)])
    if hi <= low:
        reflectionFactor = 10^10
        x = reflectionFactor - x
        newX = x
    else:
        newX = x
        x = numpy.insert(x,0,0)
        y = numpy.insert(y,0,0)
    x = numpy.array([[x1] for x1 in x])
    y = numpy.array([[y1] for y1 in y])
    lams = numpy.exp(numpy.random.rand(100, 1))
    gam = LogisticGAM(s(0,n_splines=6,spline_order=3)).gridsearch(x,y,lam=lams,progress=False)
    probs = gam.predict_proba(newX)
    odds = probs / (1 - probs)
    return(odds)
