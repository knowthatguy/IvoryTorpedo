from flask import Flask, jsonify, render_template, request, redirect, url_for
from callculation import *

app = Flask(__name__, static_url_path='/static')

@app.route("/charts", methods=["POST"])
def charts():
    r,h = map(float, (request.form.get("R"), request.form.get("H")))
    g = graphs.HGraph()
    return render_template("graph.html", fb = g.spline(r,h))

@app.route("/")
def main_form():
    return render_template("index.html")

if __name__ == '__main__':
    app.debug = True
    app.run(host='0.0.0.0')
