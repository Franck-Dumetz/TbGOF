with open("treatments.csv", "r") as f:
    text = f.read()

text = text.replace(",treated", ",__TEMP__")
text = text.replace(",untreated", ",treated")
text = text.replace(",__TEMP__", ",untreated")

with open("treatments.csv", "w") as f:
    f.write(text)
