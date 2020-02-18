/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils.taxTree;

/**
 *
 * @author bagci
 */
class TaxEdge {
    private TaxNode target;

    public TaxEdge() {
    }

    public TaxEdge(TaxNode target) {
        this.target = target;
    }

    public void setTarget(TaxNode target) {
        this.target = target;
    }

    public TaxNode getTarget() {
        return target;
    }
    
}
