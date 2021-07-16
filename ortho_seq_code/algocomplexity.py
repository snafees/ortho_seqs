d = 1
fx1d = None
fx1p = None

while d != 0:
    d = input("Enter dimensions: ")
    try:
        d = int(d)
    except:
        break
    if d == 0:
        break
    s = int(input("Enter sites: "))
    p = int(input("Enter pop size: "))

    fxd = d * p * (d * (s + 1) + 2 * s + 1)
    fxp = d * p * (d * (s + 1) + 3 * s + d * d)

    print("\nRegression will require", str(fxd), "computations.")
    if fx1d is not None:
        print(
            "Difference of",
            str(fxd - fx1d),
            "computations (" + str(round(100 * fxd / fx1d, 2)) + "% of previous value)",
        )
    fx1d = d * p * (d * (s + 1) + 2 * s + 1)
    print("\nProjection will require", str(fxp), "computations.")
    if fx1p is not None:
        print(
            "Difference of",
            str(fxp - fx1p),
            "computations (" + str(round(100 * fxp / fx1p, 2)) + "% of previous value)",
        )
    fx1p = d * p * (d * (s + 1) + 3 * s + d * d)
    print("----------------------------------------------------")
